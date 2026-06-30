# rule get_geno_trait_overlap:
#     input: 
#         genoIDs = config["inputBfile"] + ".fam",
#         traitData = "results/{trait}/{trait}_data.tsv"
#     output:
#         keepFile = "results/{trait}/IDs_keep.tsv"
#     conda:
#         "../envs/allEnv.yml"
#     shell:
#         """
#         sed '1d' {input.traitData} | awk '{{print $1}}' | fgrep -w -f - {input.genoIDs} | awk '{{print $1"\t"$2}}' > {output.keepFile}
#         """

# rule plink_filter_overlap:
#     input: 
#         keepFile = "results/{trait}/IDs_keep.tsv"
#     output:
#         bed=temp("results/{trait}/{trait}_final_geno.bed"),
#         bim=temp("results/{trait}/{trait}_final_geno.bim"),
#         fam=temp("results/{trait}/{trait}_final_geno.fam"),
#         prunedBed=temp("results/{trait}/{trait}_final_geno_pruned.bed"),
#         prunedBim=temp("results/{trait}/{trait}_final_geno_pruned.bim"),
#         prunedFam=temp("results/{trait}/{trait}_final_geno_pruned.fam")
#     params:
#         geno = config["geno"],
#         mind = config["mind"],
#         MAF = config["MAF"],
#         inputBfile = config["inputBfile"]
#     conda:
#         "../envs/allEnv.yml"
#     shell:
#         r"""
#         plink -bfile {params.inputBfile} --make-bed --keep {input.keepFile} --geno {params.geno} --mind {params.mind} --maf {params.MAF} --out results/{wildcards.trait}/{wildcards.trait}_final_geno
#         plink -bfile results/{wildcards.trait}/{wildcards.trait}_final_geno --indep-pairwise 50 5 0.2 -out results/{wildcards.trait}/{wildcards.trait}_final_geno_prune
#         plink -bfile results/{wildcards.trait}/{wildcards.trait}_final_geno --extract results/{wildcards.trait}/{wildcards.trait}_final_geno_prune.prune.in --make-bed --out results/{wildcards.trait}/{wildcards.trait}_final_geno_pruned
#         """

rule plink_filter_overlap:
    output:
        bed="results/final_geno.bed",
        bim="results/final_geno.bim",
        fam="results/final_geno.fam",
        prunedBed="results/final_geno_pruned.bed",
        prunedBim="results/final_geno_pruned.bim",
        prunedFam="results/final_geno_pruned.fam"
    params:
        geno = config["geno"],
        mind = config["mind"],
        MAF = config["MAF"],
        inputBfile = config["inputBfile"]
    conda:
        "../envs/allEnv.yml"
    shell:
        r"""
        plink -bfile {params.inputBfile} --make-bed --geno {params.geno} --mind {params.mind} --maf {params.MAF} --out results/final_geno
        plink -bfile results/final_geno --indep-pairwise 50 5 0.2 -out results/final_geno_prune
        plink -bfile results/final_geno --extract results/final_geno_prune.prune.in --make-bed --out results/final_geno_pruned
        """

# rule order_traits:
#     input:
#         genoIDs = "results/{trait}/{trait}_final_geno.fam",
#         traitData = "results/{trait}/{trait}_data.tsv"
#     output:
#         finalDat = "results/{trait}/{trait}_genoOverlap.tsv"
#     conda:
#         "../envs/allEnv.yml"
#     script:
#         "../scripts/orderTrait_byGeno.R"

rule order_traits:
    input:
        genoIDs = "results/final_geno.fam",
        traitData = "results/{trait}/{trait}_data.tsv"
    output:
        finalDat = "results/{trait}/{trait}_genoOverlap.tsv"
    conda:
        "../envs/allEnv.yml"
    script:
        "../scripts/orderTrait_byGeno.R"

# rule get_GRM:
#     input: 
#         bed="results/{trait}/{trait}_final_geno.bed",
#         bim="results/{trait}/{trait}_final_geno.bim",
#         fam="results/{trait}/{trait}_final_geno.fam"
#     output:
#         GRM = "results/{trait}/GEMMA_GRM.sXX.txt"
#     conda:
#         "../envs/allEnv.yml"
#     shell:
#         """
#         awk -F ' ' '{{if($6 == -9) {{$6 = 1}}; print $0}}' results/{wildcards.trait}/{wildcards.trait}_final_geno.fam > results/{wildcards.trait}/er3 && mv results/{wildcards.trait}/er3 results/{wildcards.trait}/{wildcards.trait}_final_geno.fam
#         gemma -bfile results/{wildcards.trait}/{wildcards.trait}_final_geno -gk 2 -o GEMMA_GRM -outdir results/{wildcards.trait}
#         """

# rule get_LOCO_GRMs:
#     input:
#         finalBim = "results/{trait}/{trait}_final_geno_pruned.bim",
#         finalBed = "results/{trait}/{trait}_final_geno_pruned.bed",
#         finalFam = "results/{trait}/{trait}_final_geno_pruned.fam"
#     output:
#         chrExclBim = temp("results/{trait}/{trait}_final_geno_nochr{chromosome}.bim"),
#         chrExclBed = temp("results/{trait}/{trait}_final_geno_nochr{chromosome}.bed"),
#         chrExclFam = temp("results/{trait}/{trait}_final_geno_nochr{chromosome}.fam"),
#         LOCOgrms = "results/{trait}/nochr{chromosome}_GRM.sXX.txt"
#     conda:
#         "../envs/allEnv.yml"
#     shell:
#         r"""
#         set -euo pipefail
#         plink --bfile results/{wildcards.trait}/{wildcards.trait}_final_geno_pruned --not-chr {wildcards.chromosome} --make-bed --out results/{wildcards.trait}/{wildcards.trait}_final_geno_nochr{wildcards.chromosome}
#         awk -F ' ' '{{if($6 == -9) {{$6 = 1}}; print $0}}' results/{wildcards.trait}/{wildcards.trait}_final_geno_nochr{wildcards.chromosome}.fam > results/{wildcards.trait}/{wildcards.trait}_final_geno_nochr{wildcards.chromosome}_tmp && mv results/{wildcards.trait}/{wildcards.trait}_final_geno_nochr{wildcards.chromosome}_tmp results/{wildcards.trait}/{wildcards.trait}_final_geno_nochr{wildcards.chromosome}.fam
#         gemma -bfile results/{wildcards.trait}/{wildcards.trait}_final_geno_nochr{wildcards.chromosome} -gk 2 -o nochr{wildcards.chromosome}_GRM -outdir results/{wildcards.trait}
#         """

rule get_LOCO_GRMs:
    input:
        finalBim = "results/final_geno_pruned.bim",
        finalBed = "results/final_geno_pruned.bed",
        finalFam = "results/final_geno_pruned.fam"
    output:
        chrExclBim = temp("results/final_geno_nochr{chromosome}.bim"),
        chrExclBed = temp("results/final_geno_nochr{chromosome}.bed"),
        chrExclFam = temp("results/final_geno_nochr{chromosome}.fam"),
        LOCOgrms = "results/nochr{chromosome}_GRM.sXX.txt"
    conda:
        "../envs/allEnv.yml"
    shell:
        r"""
        set -euo pipefail
        plink --bfile results/final_geno_pruned --not-chr {wildcards.chromosome} --make-bed --out results/final_geno_nochr{wildcards.chromosome}
        awk -F ' ' '{{if($6 == -9) {{$6 = 1}}; print $0}}' results/final_geno_nochr{wildcards.chromosome}.fam > results/final_geno_nochr{wildcards.chromosome}_tmp && mv results/final_geno_nochr{wildcards.chromosome}_tmp results/final_geno_nochr{wildcards.chromosome}.fam
        gemma -bfile results/final_geno_nochr{wildcards.chromosome} -gk 2 -o nochr{wildcards.chromosome}_GRM -outdir results
        """

# rule get_geno_PCs:
#     group: "get_PCs"
#     input:
#         bed="results/{trait}/{trait}_final_geno_pruned.bed",
#         bim="results/{trait}/{trait}_final_geno_pruned.bim",
#         fam="results/{trait}/{trait}_final_geno_pruned.fam"
#     output:
#         eigenvec="results/{trait}/{trait}_final_geno_PCA.eigenvec",
#         eigenval="results/{trait}/{trait}_final_geno_PCA.eigenval",
#         log="results/{trait}/{trait}_final_geno_PCA.log",
#         nosex="results/{trait}/{trait}_final_geno_PCA.nosex"
#     conda:
#         "../envs/allEnv.yml"
#     shell:
#         "plink --bfile results/{wildcards.trait}/{wildcards.trait}_final_geno_pruned --pca --out results/{wildcards.trait}/{wildcards.trait}_final_geno_PCA"  

rule get_geno_PCs:
    group: "get_PCs"
    input:
        bed="results/final_geno_pruned.bed",
        bim="results/final_geno_pruned.bim",
        fam="results/final_geno_pruned.fam"
    output:
        eigenvec="results/final_geno_PCA.eigenvec",
        eigenval="results/final_geno_PCA.eigenval",
        log="results/final_geno_PCA.log",
        nosex="results/final_geno_PCA.nosex"
    conda:
        "../envs/allEnv.yml"
    shell:
        "plink --bfile results/final_geno_pruned --pca --out results/final_geno_PCA"  

# rule get_final_cov:
#     group: "get_cov"
#     input: 
#         genoPCs = "results/{trait}/{trait}_final_geno_PCA.eigenvec",
#         THGPcov = config["THGPmetadata"]
#     output:
#         fullCov = "results/{trait}/{trait}_GEMMA_cov_W.tsv"
#     params:
#         genoPCnum = config["GENO_PCnum"]
#     conda:
#         "../envs/allEnv.yml"
#     script:
#         "../scripts/getFinalCovs.R"

rule get_final_cov:
    group: "get_cov"
    input: 
        genoPCs = "results/final_geno_PCA.eigenvec",
        genoIDs = "results/final_geno.fam",
        THGPcov = config["THGPmetadata"]
    output:
        fullCov = "results/final_GEMMA_cov_W.tsv"
    params:
        genoPCnum = config["GENO_PCnum"]
    conda:
        "../envs/allEnv.yml"
    script:
        "../scripts/getFinalCovs.R"

# rule get_SNP_lists:
#     input:
#         finalSNPs = "results/{trait}/{trait}_final_geno.bim"
#     output:
#         SNPlist = "results/{trait}/SNPlist_chr{chromosome}.txt"
#     shell:
#         r"""
#         chr={wildcards.chromosome}
#         awk -v chr="$chr" '$1 == chr {{ print $2 }}' {input.finalSNPs} > {output.SNPlist}
#         """

rule get_SNP_lists:
    input:
        finalSNPs = "results/final_geno.bim"
    output:
        SNPlist = "results/SNPlist_chr{chromosome}.txt"
    shell:
        r"""
        chr={wildcards.chromosome}
        awk -v chr="$chr" '$1 == chr {{ print $2 }}' {input.finalSNPs} > {output.SNPlist}
        """