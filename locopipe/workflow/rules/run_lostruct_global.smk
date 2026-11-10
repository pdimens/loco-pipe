# This rule splits chromosome-level beagle formatted genotype likelihoods into windows with fixed number of SNPs
rule split_beagle_global:
    input:
        beagle="angsd/snp_calling_global/{chr}.beagle.gz",
    output:
        done=touch("lostruct/global/split_beagle/{chr}.done"),
    params:
        outdir="lostruct/global/split_beagle",
        snp_window_size=config["lostruct"]["snp_window_size"],
    threads: 8
    shell:
        '''
        mkdir -p {params.outdir}
        rm -f {params.outdir}/{wildcards.chr}.*
        gunzip -c {input.beagle} > {params.outdir}/{wildcards.chr}.beagle
        head -n 1 {params.outdir}/{wildcards.chr}.beagle > {params.outdir}/{wildcards.chr}.beagle.header
        echo "Splitting "{wildcards.chr}
        ## Split beagle files into smaller windows, each containing a header and the desired number of SNPs
        tail -n +2 {params.outdir}/{wildcards.chr}.beagle | split -a10 -dl {params.snp_window_size} --additional-suffix=.beagle --filter='bash -c "{{ cat {params.outdir}/{wildcards.chr}.beagle.header; cat; }} | gzip > $FILE.gz"' - {params.outdir}/{wildcards.chr}.w
        rm -f {params.outdir}/{wildcards.chr}.beagle
        rm -f {params.outdir}/{wildcards.chr}.beagle.header
        '''

# This rule runs PCAngsd at each window
rule pcangsd_in_windows_global:
    input:
        done="lostruct/global/split_beagle/{chr}.done",
    output:
        done=touch("lostruct/global/run_pcangsd_in_windows/{chr}.done"),
    params:
        indir="lostruct/global/split_beagle",
        outdir="lostruct/global/run_pcangsd_in_windows",
        snp_window_size=config["lostruct"]["snp_window_size"],
        ## note that this will make the number of SNPs in each window different. if we want to keep the number of SNPs consistent, we'll need to stick with the minmaf filter used when generating the beagle file
        minmaf=MINMAF,
    threads: THREADS
    log: "lostruct/global/run_pcangsd_in_windows/{chr}.log"
    conda:
        "../envs/pcangsd.yaml" 
    shell:
        '''
        mkdir -p {params.outdir}
        rm -f {params.outdir}/{wildcards.chr}.*
        cd  {params.indir}
        for FILE in `ls {wildcards.chr}.w*.beagle.gz`; do
        N_SNPS=`zcat $FILE | wc -l`
        ## skip the last window in a chromosome if it has fewer than half of the SNPs in other windows
        if [ $N_SNPS -gt $(({params.snp_window_size}/2)) ]; then
        pcangsd \
        --beagle $FILE \
        --maf {params.minmaf} \
        --threads {threads} \
        --out {params.outdir}/${{FILE%%.beagle.gz}} \
        &>> {log}
        fi
        done
        '''

# This rule summarizes PCA results at each window as required by the lostruct package
rule summarize_pcangsd_for_lostruct_global:
    input:
        done="lostruct/global/run_pcangsd_in_windows/{chr}.done",
    output:
        pca_summary="lostruct/global/summarize_pcangsd_for_lostruct/{chr}.pca_summary.tsv",
        snp_position="lostruct/global/summarize_pcangsd_for_lostruct/{chr}.snp_position.tsv",
        done=touch("lostruct/global/summarize_pcangsd_for_lostruct/{chr}.done"),
    params:
        outdir="lostruct/global/summarize_pcangsd_for_lostruct",
        cov_dir="lostruct/global/run_pcangsd_in_windows",
        beagle_dir="lostruct/global/split_beagle",
        rscript=config["scriptdir"] + "/summarize_pcangsd_for_lostruct.R",
        pc=config["lostruct"]["pc"],
    threads: 8
    log: "lostruct/global/summarize_pcangsd_for_lostruct/{chr}.log"
    conda:
        "lostruct_lcpipe"  ## lostruct isn't available on conda. we'll need to set up a conda environment with lostruct beforehand
    shell:
        '''
        mkdir -p {params.outdir}
        rm -f {params.outdir}/{wildcards.chr}.*
        cd  {params.cov_dir}
        for FILE in `ls {wildcards.chr}.w*.cov`; do
        Rscript --vanilla {params.rscript} ${{FILE%%.cov}} {wildcards.chr} {params.pc} {params.cov_dir} {params.beagle_dir} {params.outdir} &>> {log}
        done
        '''

# This rule runs the lostruct package using PCA summaries as input and outputs a distance matrix, which describes the dissimilarity between windows
rule lostruct_global:
    input:
        done=expand("lostruct/global/summarize_pcangsd_for_lostruct/{chr}.done", chr=CHRS),
    output:
        dist="lostruct/global/run_lostruct/combined.dist.tsv",
        mds="lostruct/global/run_lostruct/combined.mds.tsv",
        proportion_variance="lostruct/global/run_lostruct/combined.proportion_variance.txt",
        mds_by_chr=expand("lostruct/global/run_lostruct/{chr}.mds.tsv", chr=CHRS),
        proportion_variance_by_chr=expand("lostruct/global/run_lostruct/{chr}.proportion_variance.txt", chr=CHRS),
        done=touch("lostruct/global/run_lostruct/run_lostruct.done"),
    params:
        rscript=config["scriptdir"] + "/run_lostruct.R",
        chr_table_path="workflow/" + config["chr_table"],
        pc=config["lostruct"]["pc"],
        k=config["lostruct"]["k"],
        indir="lostruct/global/summarize_pcangsd_for_lostruct",
        outdir="lostruct/global/run_lostruct",
    threads: THREADS
    log: "lostruct/global/run_lostruct/run_lostruct.log"
    conda:
        "lostruct_lcpipe"  ## lostruct isn't available on conda. we'll need to set up a conda environment with lostruct beforehand
    shell:
        '''
        mkdir -p {params.outdir}
        rm -f {params.outdir}/*
        Rscript --vanilla {params.rscript} {params.chr_table_path} {params.pc} {params.k} {params.indir} {params.outdir} {threads} &> {log}
        '''

# This rule plots the top MDS axes and saves the outlier windows in text files
rule plot_lostruct_mds_global:
    input:
        done="lostruct/global/run_lostruct/run_lostruct.done",
    output:
        combined_table="lostruct/global/plot_lostruct_mds/combined.mds.csv",
        combined_plot="figures/lostruct/global/combined.mds.png",
        separated_table="lostruct/global/plot_lostruct_mds/separated.mds.csv",
        separated_plot="figures/lostruct/global/separated.mds.png",
        done=touch("figures/lostruct/global/plot_lostruct_mds.done"),
    params:
        rscript=config["scriptdir"] + "/plot_lostruct_mds.R",
        chr_table_path="workflow/" + config["chr_table"],
        indir="lostruct/global/run_lostruct",
        outdir="lostruct/global/plot_lostruct_mds",
        plot_dir="figures/lostruct/global",
        k=config["lostruct"]["k"],
        z_cutoff=config["lostruct"]["z_cutoff"],
        fig_width=config["lostruct"]["fig_width"],
        fig_height=config["lostruct"]["fig_height"],
    threads: 4
    log: "lostruct/global/plot_lostruct_mds/plot_lostruct_mds.log"
    conda:
        "lostruct_lcpipe"  ## lostruct isn't available on conda. we'll need to set up a conda environment with lostruct beforehand
    shell:
        '''
        mkdir -p {params.plot_dir}
        mkdir -p {params.outdir}
        rm -f {params.outdir}/*
        Rscript --vanilla {params.rscript} {params.chr_table_path} {params.indir} {params.outdir} {params.plot_dir} {params.k} {params.z_cutoff} {params.fig_width} {params.fig_height} &> {log}
        '''

# This rule combines outlier windows along the same MDS axis and performs a concensus PCA with PCAngsd
rule pcangsd_with_lostruct_outliers_global:
    input:
        done="figures/lostruct/global/plot_lostruct_mds.done",
    output:
        done=touch("lostruct/global/run_pcangsd_with_lostruct_outliers/run_pcangsd_with_lostruct_outliers.done"),
    params:
        indir="lostruct/global/plot_lostruct_mds",
        outdir="lostruct/global/run_pcangsd_with_lostruct_outliers",
        beagle_dir="lostruct/global/split_beagle",
        minmaf=MINMAF,
    threads: THREADS
    log: "lostruct/global/run_pcangsd_with_lostruct_outliers/run_pcangsd_with_lostruct_outliers.log"
    conda:
        "../envs/pcangsd.yaml" 
    shell:
        '''
        mkdir -p {params.outdir}
        rm -f {params.outdir}/*
        cd  {params.indir}
        ## Iterate over each MDA axis
        for FILE in `ls *.tsv`; do
        ## Assemble beagle files for each MDS axis
        PREFIX=${{FILE%%.tsv}}
        I=1
        ## Iterate over each window
        for WINDOW in `cut -f 1 $FILE`; do
        if [ $I -eq 1 ] ; then
            zcat {params.beagle_dir}/${{WINDOW}}.beagle.gz > {params.outdir}/${{PREFIX}}.beagle
        else zcat {params.beagle_dir}/${{WINDOW}}.beagle.gz | tail -n +2 >> {params.outdir}/${{PREFIX}}.beagle
        fi
        I=$(( I + 1 ))
        done
        gzip {params.outdir}/${{PREFIX}}.beagle
        ## Run PCAngsd
        pcangsd \
        --beagle {params.outdir}/${{PREFIX}}.beagle.gz \
        --maf {params.minmaf} \
        --threads {threads} \
        --out {params.outdir}/${{PREFIX}} \
        &>> {log}
        done
        '''

# This rule plots the result of the outlier windows' concensus PCA
rule plot_lostruct_outlier_pca_global:
    input:
        done="lostruct/global/run_pcangsd_with_lostruct_outliers/run_pcangsd_with_lostruct_outliers.done",
    output:
        combined_table="lostruct/global/plot_lostruct_outlier_pca/combined.pca.tsv",
        combined_plot="figures/lostruct/global/combined.pca.png",
        separated_table="lostruct/global/plot_lostruct_outlier_pca/separated.pca.tsv",
        separated_plot="figures/lostruct/global/separated.pca.pdf",
        done=touch("figures/lostruct/global/plot_lostruct_outlier_pca.done"),
    params:
        outdir="lostruct/global/plot_lostruct_outlier_pca",
        rscript=config["scriptdir"] + "/plot_lostruct_outlier_pca.R",
        cov_dir="lostruct/global/run_pcangsd_with_lostruct_outliers",
        plot_dir="figures/lostruct/global",
        sample_table_path="workflow/" + config["sample_table"],
        color_by=POP_L1_COL,
        chr_table_path="workflow/" + config["chr_table"],
        k=config["lostruct"]["k"],
    threads: 4
    log: "lostruct/global/plot_lostruct_outlier_pca/plot_lostruct_outlier_pca.log"
    conda:
        "lostruct_lcpipe"  ## lostruct isn't available on conda. we'll need to set up a conda environment with lostruct beforehand
    shell:
        '''
        mkdir -p {params.outdir}
        Rscript --vanilla {params.rscript} {params.cov_dir} {params.plot_dir} {params.sample_table_path} {params.color_by} {params.chr_table_path} {params.k} &> {log}
        '''
