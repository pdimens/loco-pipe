# This rule takes in population-specific bamlists and then estimate allele frequencies at each SNP location without subsetting. It will also output population-specific genotype likelihood in Beagle format and site allele frequency likelihoods.
rule get_maf:
    input:
        ref= REFERENCE,
        fai= REFERENCE+".fai",
        bamlist="workflow/{population}_" + BAMLIST,
        site_list="angsd/snp_calling_global/{chr}.snp_list"
    output:
        expand("angsd/get_maf/{{population}}.{{chr}}.{ext}",
               ext=[ "arg", "beagle.gz", "covMat", "ibs.gz", "ibsMat", "mafs.gz", "pos.gz", "saf.gz", "saf.idx", "saf.pos.gz" ]),
        touch("angsd/get_maf/{population}.{chr}.done")
    threads:
        THREADS
    params:
        outdir="angsd/get_maf",
        gl_model=config["gl_model"],
        minq=MINQ,
        minmapq=MINMAPQ,
        extra=config["get_maf"]["extra"],
        minind=MININD,
        mindepthind=MINDEPTHIND,
    log: "angsd/get_maf/{population}.{chr}.log"
    conda: "../envs/angsd.yaml"
    shell:
        '''
        mkdir -p {params.outdir}
        angsd \
            -bam {input.bamlist} -anc {input.ref} \
            -P {threads} \
            -out {params.outdir}/{wildcards.population}.{wildcards.chr} \
            -sites {input.site_list} -r {wildcards.chr} \
            -GL {params.gl_model} -doGlf 2 -doMaf 1 -doSaf 1 -doMajorMinor 3 \
            -minInd {params.minind} -setMinDepthInd {params.mindepthind} \
            -doCounts 1 -dumpCounts 1 \
            -doIBS 1 -makematrix 1 -doCov 1 \
            -minQ {params.minq} -minMapQ {params.minmapq} \
            -remove_bads 1 -only_proper_pairs 1 \
            {params.extra} 2> {log}
        '''
