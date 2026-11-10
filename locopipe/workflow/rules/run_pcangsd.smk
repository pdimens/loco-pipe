# This rule uses PCAngsd software to conduct PCA analysis on a global level. 
# More details about PCAngsd could be found here: http://www.popgen.dk/software/index.php/PCAngsd.
rule run_pcangsd_pca_global:
    input: "angsd/snp_calling_global/{file}.beagle.gz",
    output: 
        cov = "pcangsd/global/{file}.cov",
        done = touch("pcangsd/global/{file}.done"),
    params:
        outdir = "pcangsd/global",
        minmaf=MINMAF,
    threads: THREADS
    log: "pcangsd/global/{file}.log"
    conda:
        "../envs/pcangsd.yaml" 
    shell:
        '''
        conda env list > {log}
        mkdir -p {params.outdir}
        ## run pcangsd for PCA only (other options need to be added)
        pcangsd \
        --beagle {input} \
        --snp-weights \
        --sites-save \
        --maf {params.minmaf} \
        --threads {threads} \
        --out {params.outdir}/{wildcards.file} \
        &>> {log}
        '''

# This rule uses PCAngsd software to conduct PCA analysis on a local level. 
# More details about PCAngsd could be found here: http://www.popgen.dk/software/index.php/PCAngsd.      
rule run_pcangsd_pca_local:
    input: "angsd/get_maf/{population}.{file}.beagle.gz",
    output: 
        cov = "pcangsd/local/{population}.{file}.cov",
        done = touch("pcangsd/local/{population}.{file}.done"),
    params:
        outdir = "pcangsd/local",
        minmaf = MINMAF,
    threads: THREADS
    log: "pcangsd/local/{population}.{file}.log"
    conda:
        "../envs/pcangsd.yaml" 
    shell:
        '''
        mkdir -p {params.outdir}
        ## run pcangsd for PCA only (other options need to be added)
        pcangsd \
        --beagle {input} \
        --snp-weights \
        --sites-save \
        --maf {params.minmaf} \
        --threads {threads} \
        --out {params.outdir}/{wildcards.population}.{wildcards.file} \
        &> {log}
        '''

# This rule plots the PCA analysis conducted in the rule,run_pcangsd_pca_global, above. 
rule plot_pcangsd_pca_global:
    input:
        cov = "pcangsd/global/{file}.cov",
        done = "pcangsd/global/{file}.done",
    output:
        table = "pcangsd/global/{file}.tsv",
        plot = "figures/pcangsd/global/{file}.png",
        done = touch("figures/pcangsd/global/{file}.done"),
    params:
        outdir = "figures/pcangsd/global",
        sample_table_path = "workflow/" + config["sample_table"],
        color_by = POP_L1_COLNAME,
        rscript = config["scriptdir"] + "/plot_pcangsd_pca.R",
    threads: 1
    log: "pcangsd/global/{file}.plot.log"
    conda:
        "../envs/r.yaml" 
    shell:
        '''
        mkdir -p {params.outdir}
        Rscript --vanilla {params.rscript} {input.cov} {output.table} {output.plot} {params.sample_table_path} {params.color_by} &> {log}
        '''

# This rule plots the PCA analysis conducted in the rule,run_pcangsd_pca_local, above. 
rule plot_pcangsd_pca_local:
    input:
        cov = "pcangsd/local/{population}.{file}.cov",
        done = "pcangsd/local/{population}.{file}.done",
    output:
        table = "pcangsd/local/{population}.{file}.tsv",
        plot = "figures/pcangsd/local/{population}.{file}.png",
        done = touch("figures/pcangsd/local/{population}.{file}.done"),
    params:
        outdir = "figures/pcangsd/local",
        sample_table_path = "workflow/" + config["sample_table"],
        color_by = POP_L1_COLNAME,
        pop_col = POP_L1_COLNAME,
        rscript = config["scriptdir"] + "/plot_pcangsd_pca.R",
    threads: 1
    log: "pcangsd/local/{population}.{file}.plot.log"
    conda:
        "../envs/r.yaml" 
    shell:
        '''
        mkdir -p {params.outdir}
        Rscript --vanilla {params.rscript} {input.cov} {output.table} {output.plot} {params.sample_table_path} {params.color_by} {params.pop_col} {wildcards.population} &> {log}
        '''
