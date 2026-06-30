# rule get_GEMMA_files:
#     group: "get_plink"
#     input: 
#         finalDat = "results/{trait}/{trait}_genoOverlap.tsv",
#         finalFam = "results/{trait}/{trait}_final_geno.fam",
#         finalBed = "results/{trait}/{trait}_final_geno.bed",
#         finalBim = "results/{trait}/{trait}_final_geno.bim"
#     output:
#         traitFam = "results/{trait}/{trait}.fam",
#         traitBed = temp("results/{trait}/{trait}.bed"),
#         traitBim = temp("results/{trait}/{trait}.bim")
#     shell:
#         """
#         cat {input.finalDat} | paste {input.finalFam} - -d' ' | cut -f1-5,7 -d ' ' > {output.traitFam}
#         cp {input.finalBed} {output.traitBed}
#         cp {input.finalBim} {output.traitBim}
#         """

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

# rule run_ME:
#     input:
#         fullCov = "results/{trait}/{trait}_GEMMA_cov_W.tsv",
#         plinkFiles = ["results/{trait}/{trait}.bim",
#                       "results/{trait}/{trait}.bed",
#                       "results/{trait}/{trait}.fam"],
#         LOCOgrms = "results/{trait}/nochr{chromosome}_GRM.sXX.txt",
#         SNPlist = "results/{trait}/SNPlist_chr{chromosome}.txt"
#     output:
#         traitRes = "results/{trait}/{trait}_chr{chromosome}_LMM.assoc.txt"
#     conda:
#         "../envs/allEnv.yml"
#     shell:
#         r"""
#         set -euo pipefail
#         thisChr={wildcards.chromosome}

#         gemma \
#             -bfile results/{wildcards.trait}/{wildcards.trait} \
#             -k results/{wildcards.trait}/nochr${{thisChr}}_GRM.sXX.txt \
#             -c {input.fullCov} \
#             -snps {input.SNPlist} \
#             -miss 1.0 -lmm 4 \
#             -outdir results/{wildcards.trait} \
#             -o {wildcards.trait}_chr${{thisChr}}_LMM
#         """

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

# rule cleanup:
#     input:
#         bfiles_orig = expand("results/{trait}/{trait}_final_geno.{ext}", trait = traits, ext = ["bim", "bed", "ped", "map", "fam"]),
#         bfiles = expand("results/{trait}/{trait}.{ext}", trait = traits, ext = ["bim", "bed", "fam"])
#     shell:
#         """
#         rm -f {input.bfiles_orig}
#         rm -f {input.bfiles}
#         """