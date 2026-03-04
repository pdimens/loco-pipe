import click
from importlib import resources
from itertools import chain
import os
from pathlib import Path
import shutil
import subprocess
import sys
import yaml
from locopipe.workflow import CONFIG
from locopipe.validation import BAMfile, valid_dir

@click.help_option('--help', hidden = True)
@click.group()
def cli():
    """Convenience features to make loco-pipe friendlier"""

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False})
@click.option('-o', "--output", default = ".", show_default = True, type = str, required = False, help = "name of directory to create")
@click.option("-r", "--reference", required = True, type = click.Path(exists=True, dir_okay=False, readable=True), help = "reference fasta file which BAMs were aligned to")
@click.argument("samples", nargs = -1, required = True, type = BAMfile())
@click.help_option('--help', hidden = True)
def init(reference, output, samples):
    """
    Initialize a configured loco-pipe project directory

    Create a directory (--output) with all the components loco-pipe needs to run. Must include
    a `--reference` file. Provide the input BAM files and/or directories at the end of the command
    as individual files/folders, using shell wildcards (e.g. `data/acronotus*.fq`), or both.
    """
    #dir = os.path.join(output, "workflow")
    os.makedirs(os.path.join(output, "docs"), exist_ok=True)
    notices = []
    shutil.copytree(
        resources.files('locopipe.workflow').as_posix(), "workflow",
        dirs_exist_ok=True,
        ignore = shutil.ignore_patterns('*__init__.py', '*__pycache__')    
    )
    profile = {
        "directory": os.path.abspath(output),
        "rerun-incomplete": True,
        "show-failed-logs": True,
        "rerun-triggers": ["mtime", "params"],
        "quiet": "reason",
        "scheduler": "greedy",
        "software-deployment-method": "conda",
        "conda-cleanup-pkgs": "cache"
    }
    with open(os.path.join("workflow", 'config.yaml'), "w", encoding="utf-8") as sm_config:
        yaml.dump(profile, sm_config, sort_keys=False)

    with open(os.path.join(output, "locopipe.yaml"), 'w') as cnfg:
        cnf = CONFIG.format(
            basedir = os.path.abspath(output),
            reference = os.path.abspath(reference),
            scriptdir = os.path.abspath("workflow/scripts"),
            groupcol = "population",
            sampletable = "samples.tsv",
            chromtable = "contigs.tsv"
        )
        cnfg.write(cnf + "\n")
    
    with open(reference, 'r') as fa, open(os.path.join(output, "docs", "contigs.tsv"), 'w') as tb:
        id = 0
        for i in fa:
            if i.startswith(">"):
                id += 1
                _name = i.lstrip(">").split(" ")[0].strip()
                tb.write(f"{_name}\t{id}\n")
    
    if id >= 100:
        notices.append(
            "- Number of included chromosomes is or exceeds 100, which is not recommended. "
            f"The chromosome table ({os.path.join(output, 'contigs.tsv')}) should be curated to include <100."
        )
    with open(os.path.join(output, "docs", "samples.tsv"), 'w') as tb:
        tb.write("sample_name\tbam\tpopulation\n")
        for i in list(chain.from_iterable(samples)):
            _name = Path(i).stem
            tb.write(f"{_name}\t{i}\tgroupname\n")
    
    notices.append(
        f"- Samples have all been assigned to a single group and need to have distinct classifications in {os.path.join(output, 'samples.tsv')}."
    )

    print("Notices:\n" + "\n".join(notices), file = sys.stderr)

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False})
@click.option('-@', '--threads', default = 8, show_default = True, type = click.IntRange(1,999, clamp = True), help = 'Number of threads to use')
@click.argument("directory", default = ".", type = click.Path(exists=True, file_okay=False, readable=True))
@click.help_option('--help', hidden = True)
def pipe(directory, threads):
    """
    Launch the loco-pipe snakemake workflow
    
    The target `directory` is expected to have been created with `loco init`.
    Make sure you have addressed any notices provided by the init command prior to launching
    loco-pipe. Will print snakemake output to the terminal.
    """
    prjdir = os.path.join(os.path.relpath("workflow"), "pipelines", "loco-pipe.smk")
    #absent = valid_dir(directory)
    #if absent:
    #    print("Error: target directory is missing necessary loco-pipe files:\n" + "\n".join(absent), file = sys.stderr)
    #    sys.exit(1)
    cmd = ['snakemake', '--cores', f"{threads}"]
    cmd += ["--snakefile", prjdir]
    cmd += ["--configfile", os.path.join(os.path.relpath(directory), "locopipe.yaml")]
    cmd += ["--profile", os.path.relpath("workflow")]
    print("Snakemake command:\n" + " ".join(cmd))
    try:
        subprocess.run(cmd)
    except KeyboardInterrupt:
        pass

cli.add_command(pipe)
cli.add_command(init)