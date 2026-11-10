# This snakemake file runs Ohana to conduct admixture analysis to characterize population structure. 
# More detailed description of Ohana could be found here: https://github.com/jade-cheng/ohana/blob/master/README.md.

# This rule converts the subsetted beagle formatted genotype likelihood file to an .lgm format file 
# which is needed for Ohana.
rule get_ohana_input_global:
    input: "angsd/snp_calling_global/{file}.beagle"
    output: "ohana/global/{file}.lgm"
    params: outdir = "ohana/global"
    conda: "../envs/ohana.yaml"
    threads: 1
    shell:
        '''
        mkdir -p {params.outdir}
        convert bgl2lgm {input} {output}
        '''
        
# After the file conversion from the get_ohana_input_global rule, this rule applies "qpas"
# method in Ohana with different hyperparameters for k. You could adjust the range for k by 
# changing the "min_k" and "max_k" in the config file. The final outputs are q_matrix and f_matrix. 
# The former records the admixture proportions for each sample, and the latter records the allele 
# frequencies in each ancestral population.
# A well-explained example could be found here: https://github.com/jade-cheng/ohana/blob/master/README.md.
rule ohana_global:
     input: "ohana/global/{file}.lgm"
     output: 
         q_mat = "ohana/global/{file}.k{k}.q.matrix",
         f_mat = "ohana/global/{file}.k{k}.f.matrix",
         done = touch("ohana/global/{file}.k{k}.done"),
     threads: 4
     log: "ohana/global/{file}.k{k}.log"
     conda: "../envs/ohana.yaml"
     shell:
         '''
         qpas {input} \
         -k {wildcards.k} \
         -qo {output.q_mat} \
         -fo {output.f_mat} \
         -mi 50 &> {log}
         '''
# The next two rules inherit the same logic from get_ohana_input_global and run_ohana_global rules,
# but these rules uses data specific to one population.
rule get_ohana_input_local:
    input: "angsd/get_maf/{population}.{file}.beagle"
    output: "ohana/local/{population}.{file}.lgm"
    params: outdir = "ohana/local"
    conda: "../envs/ohana.yaml"
    threads: 1
    shell:
        '''
        mkdir -p {params.outdir}
        convert bgl2lgm {input} {output}
        '''

rule ohana_local:
     input: "ohana/local/{population}.{file}.lgm"
     output: 
         q_mat = "ohana/local/{population}.{file}.k{k}.q.matrix",
         f_mat = "ohana/local/{population}.{file}.k{k}.f.matrix",
         done = touch("ohana/local/{population}.{file}.k{k}.done"),
     threads: 4
     log: "ohana/local/{population}.{file}.k{k}.log"
     conda: "../envs/ohana.yaml"
     shell:
         '''
         qpas {input} \
         -k {wildcards.k} \
         -qo {output.q_mat} \
         -fo {output.f_mat} \
         -mi 50 &> {log}
         '''

rule plot_ohana_admixture_global:
    input:
        done = expand("ohana/global/{{file}}.k{k}.done", k=list(range(config["ohana_global"]["min_k"], config["ohana_global"]["max_k"] + 1))),
    output:
        table = "ohana/global/{file}.tsv",
        plot = "figures/ohana/global/{file}.png",
        done = touch("figures/ohana/global/{file}.done"),
    params:
        indir = "ohana/global",
        outdir = "figures/ohana/global",
        sample_table_path = "workflow/" + config["sample_table"],
        group_by = POP_L1_COL,
        min_k = config["ohana_global"]["min_k"],
        max_k = config["ohana_global"]["max_k"],
        rscript = config["scriptdir"] + "/plot_ohana_admixture.R",
    threads: 1
    log: "ohana/global/{file}.plot.log"
    conda:
        "../envs/r.yaml" 
    shell:
        '''
        mkdir -p {params.outdir}
        Rscript --vanilla {params.rscript} {params.indir} {wildcards.file} {output.table} {output.plot} {params.sample_table_path} {params.group_by} {params.min_k} {params.max_k} &> {log}
        '''

rule plot_ohana_admixture_local:
    input:
        done = expand("ohana/local/{{population}}.{{file}}.k{k}.done", k=list(range(config["ohana_local"]["min_k"], config["ohana_local"]["max_k"] + 1))),
    output:
        table = "ohana/local/{population}.{file}.tsv",
        plot = "figures/ohana/local/{population}.{file}.png",
        done = touch("figures/ohana/local/{population}.{file}.done"),
    params:
        indir = "ohana/local",
        outdir = "figures/ohana/local",
        sample_table_path = "workflow/" + config["sample_table"],
        group_by = POP_L1_COL,
        min_k = config["ohana_local"]["min_k"],
        max_k = config["ohana_local"]["max_k"],
        pop_col = POP_L1_COLNAME,
        rscript = config["scriptdir"] + "/plot_ohana_admixture.R",
    threads: 1
    log: "ohana/local/{population}.{file}.plot.log"
    conda:
        "../envs/r.yaml" 
    shell:
        '''
        mkdir -p {params.outdir}
        Rscript --vanilla {params.rscript} {params.indir} {wildcards.file} {output.table} {output.plot} {params.sample_table_path} {params.group_by} {params.min_k} {params.max_k} {params.pop_col} {wildcards.population} &> {log}
        '''
