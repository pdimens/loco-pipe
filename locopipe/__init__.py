import click
from importlib import resources
from itertools import chain
import os
from pathlib import Path
import shutil
import subprocess
import sys
import yaml
from locopipe.validation import BAMfile

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False})
@click.option('-s', '--simple', is_flag = True, default = False, show_default = True, help = 'Do not annotate the workflow YAML file')
@click.option('--hpc', type = click.Choice(['slurm', 'lsf', 'sge'], case_sensitive=False), help = 'Include an HPC profile for this scheduler')
@click.argument('output', type = click.Path(exists = False, resolve_path = True), required = True)
@click.argument('reference', required = True, type = click.Path(exists=True, dir_okay=False, readable=True))
@click.argument("samples", nargs = -1, required = True, type = BAMfile())
@click.help_option('--help', hidden = True)
def init(output, reference, samples, simple, hpc):
    """
    Initialize a configured loco-pipe project

    Create an OUTPUT directory with all the components loco-pipe needs to run. Provide a REFERENCE
    FASTA file that was used to create the alignment files. Provide the input BAM files and/or
    directories at the end of the command as individual files/folders, using shell wildcards (e.g. `data/acro*.bam`), or both.
    """
    checklist = []
    ## ----- Establish workflow directory --------------------
    workdir = os.path.relpath(os.path.join(output, "workflow"))
    os.makedirs(workdir, exist_ok=True)

    ## ----- Copy over loco-pipe -----------------------------
    shutil.copytree(
        resources.files('locopipe.workflow').as_posix(), f"{output}/workflow",
        dirs_exist_ok=True,
        ignore = shutil.ignore_patterns('*__init__.py', '*__pycache__', 'configs')    
    )

    ## ----- Setup and write profile.yaml -------------------
    profile = {
        **({'__use_yte__': True} if hpc else {}),
        "directory": os.path.abspath(output),
        "rerun-incomplete": True,
        "show-failed-logs": True,
        "rerun-triggers": ["mtime", "params"],
        "quiet": "reason",
        "scheduler": "greedy",
        "software-deployment-method": "conda",
        "conda-cleanup-pkgs": "cache"
    }
    if hpc:
        _hpcyaml = resources.files("locopipe.configs") / f'{hpc.lower()}.yml'
        with open(_hpcyaml) as f:
            _hpc = yaml.full_load(f)
        # merge the dicts
        profile = profile | _hpc
        checklist.append(f'• install \033[92msnakemake-executor-plugin-{hpc}\033[0m (if needed)')
        checklist.append(
            f'• provide HPC configurations to \033[92m{os.path.join(workdir, 'profile.yaml')}\033[0m'
        )

    with open(os.path.join(workdir, 'profile.yaml'), "w", encoding="utf-8") as sm_profile:
        yaml.dump(profile, sm_profile, sort_keys=False)


    ## ----- Setup and write locopipe.yaml, the config yaml file -------------- 
    with open(os.path.join(workdir, "locopipe.yaml"), 'w') as cnfg:
        with open(resources.files("locopipe.configs") / 'config.yaml', 'r') as f:
            cnf = f.read().format(
                basedir = os.path.abspath(output),
                reference = os.path.abspath(reference),
                scriptdir = os.path.abspath(f"{output}/workflow/scripts"),
                groupcol = "population",
                sampletable = "samples.tsv",
                chromtable = "contigs.tsv"
            )
        if not simple:
            cnfg.write(cnf + "\n")
        else:
            for line in cnf.splitlines():
                if line.strip() and '#' not in line:  
                  cnfg.write(line + '\n')
        checklist.append(
            f'• adjust \033[92m{os.path.join(workdir, 'locopipe.yaml')}\033[0m for your project'
        )


    ## ----- Populate contigs.tsv with numbered contig names ---------
    try:
        with open(reference, 'r') as fa:
            # test by reading 1 byte
            fa.read(1)
        fa = open(reference, 'r')
    except UnicodeDecodeError:
        import gzip
        fa = gzip.open(reference, 'rt')
    with open(os.path.join(workdir, "contigs.tsv"), 'w') as tb:
        id = 0
        for i in fa:
            if i.startswith(">"):
                id += 1
                _name = i.lstrip(">").split(" ")[0].strip()
                tb.write(f"{_name}\t{id}\n")
    fa.close()

    if id >= 80:
        checklist.append(f"• curate \033[92m{os.path.join(os.path.relpath(output), "docs", "contigs.tsv")}\033[0m so it has <80 contigs")


    ## ----- Setup samples.tsv with sample names and BAM files ------------
    with open(os.path.join(workdir, "samples.tsv"), 'w') as tb:
        tb.write("sample_name\tbam\tpopulation\n")
        for i in list(chain.from_iterable(samples)):
            _name = Path(i).stem
            tb.write(f"{_name}\t{i}\tgroupname\n")
    
    checklist.append(
        f'• assign samples to populations in \033[92m{os.path.join(os.path.relpath(output), 'docs', 'samples.tsv')}\033[0m'
    )


    ## ----- Print checklist based on above ---------------
    print("\nChecklist before running loco-pipe", file = sys.stderr)
    print("──────────────────────────────────", file = sys.stderr)
    print("\n".join(checklist), file = sys.stderr)

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False})
@click.option('-t', '--threads', default = 8, show_default = True, type = click.IntRange(1,999, clamp = True), help = 'Number of threads to use')
@click.option('-d', '--dry', is_flag = True, default = False, show_default = True, help = 'Perform a snakemake dry run')
@click.argument("directory", default = ".", type = click.Path(exists=True, file_okay=False, readable=True))
@click.help_option('--help', hidden = True)
def start(directory, threads, dry):
    """
    Launch the loco-pipe snakemake workflow
    
    The target `directory` is expected to have been created with `loco init`.
    Make sure you have addressed any notices provided by the init command prior to launching
    loco-pipe. Will print snakemake output to the terminal.
    """
    workflowdir = os.path.relpath(os.path.join(directory, 'workflow'))
    snakefile = os.path.join(workflowdir, "pipelines", "loco-pipe.smk")
    cmd = ['snakemake', '--cores', f"{threads}"]
    cmd += ["--snakefile", snakefile]
    cmd += ["--configfile", os.path.join(workflowdir, "locopipe.yaml")]
    cmd += ["--profile", os.path.join(workflowdir, 'profile.yaml')]
    if dry:
        cmd += ["-n"]
    print("Starting loco-pipe...", file = sys.stderr)
    try:
        subprocess.run(cmd)
    except KeyboardInterrupt:
        pass
