rule get_GEMMA_files:
    group: "get_plink"
    input: 
        finalDat = "results/{trait}/{trait}_genoOverlap.tsv",
        finalFam = "results/final_geno.fam",
        finalBed = "results/final_geno.bed",
        finalBim = "results/final_geno.bim"
    output:
        traitFam = "results/{trait}/{trait}.fam",
        traitBed = temp("results/{trait}/{trait}.bed"),
        traitBim = temp("results/{trait}/{trait}.bim")
    shell:
        """
        cat {input.finalDat} | paste {input.finalFam} - -d' ' | cut -f1-5,7 -d ' ' > {output.traitFam}
        cp {input.finalBed} {output.traitBed}
        cp {input.finalBim} {output.traitBim}
        """

rule run_ME:
    input:
        fullCov = "results/final_GEMMA_cov_W.tsv",
        plinkFiles = ["results/{trait}/{trait}.bim",
                      "results/{trait}/{trait}.bed",
                      "results/{trait}/{trait}.fam"],
        LOCOgrms = "results/nochr{chromosome}_GRM.sXX.txt",
        SNPlist = "results/SNPlist_chr{chromosome}.txt"
    output:
        traitRes = "results/{trait}/{trait}_chr{chromosome}_LMM.assoc.txt"
    conda:
        "../envs/allEnv.yml"
    shell:
        r"""
        set -euo pipefail
        thisChr={wildcards.chromosome}

        gemma \
            -bfile results/{wildcards.trait}/{wildcards.trait} \
            -k results/{wildcards.trait}/nochr${{thisChr}}_GRM.sXX.txt \
            -c {input.fullCov} \
            -snps {input.SNPlist} \
            -miss 1.0 -lmm 4 \
            -outdir results/{wildcards.trait} \
            -o {wildcards.trait}_chr${{thisChr}}_LMM
        """
#         """
#         rm -f {input.bfiles_orig}
#         rm -f {input.bfiles}
#         """
