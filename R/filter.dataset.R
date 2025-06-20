#' @title Quality control of genotype call files (VCF, BCF, PLINK files)
#' @description
#' This function performs quality control of large-scale datasets with multiple filters.
#' @param input.file (character) The file path to the genetic dataset. Accepts vcf or bcf files.
#' @param plink.files (logical) If TRUE, PLINK files will be used (bed, bim, fam).
#' @param bed.file (character) The file path for the bed file. Should be provided if plink.files = TRUE
#' @param bim.file (character) The file path for the bim file. Should be provided if plink.files = TRUE
#' @param fam.file (character) The file path for the fam file. Should be provided if plink.files = TRUE
#' @param remove.related (logical) If TRUE, function will calculate the kinship coefficient (to be supplied) of samples and remove those statistically likely to be related. Default is FALSE.
#' @param kinship.coefficient (optional) The kinship threshold. Based on kingrelatedness.com/manual, a kinship coefficient range of >0.354 corresponds to duplicate/MZ twins, [0.177, 0.354] to 1st-degree relationships, [0.0844, 0.177] to 2nd-degree relationships, and [0.0442, 0.0884] to 3rd-degree relationships.
#' @param geno.value (integer) The percentage threshold of individuals with a missing genotypes to be filtered.
#' @param maf.value (integer) The minor allele frequency threshold. 
#' @param mind.value (integer) The percentage threshold of individuals with missing data to be removed.
#' @param limit.LD (logical) If TRUE, function should prune the markers that are in linkage equilibrium. Default is FALSE
#' @param indep.pairwise.kb (integer) should be set if limit.LD is TRUE. This is the window size in variant count or kilobase units. 
#' @param indep.pairwise.ct (integer) should be set if limit.LD is TRUE. This is the variant count to shift the window at each end of step.
#' @param indep.pairwise.r2 (integer) should be set if limit.LD is TRUE. This is the pairwise threshold. At each step, pairs of variants in the current window with squared correlation greater than the threshold are noted and the variants are pruned from the window until no such pairs remain.
#' @examples
#' filter.dataset("test_data/test.csv", plink.files = FALSE, remove.related = TRUE, kinship.coefficient = 0.0844, geno.value = 0.5, maf.value = 0.05, mind.value = 0.1, limit.LD = FALSE)
#' filter.dataset("test_data/kgpdata.vcf", plink.files = FALSE, remove.related = FALSE, geno.value = 0.5, maf.value = 0.05, mind.value = 0.1, limit.LD = TRUE, indep.pairwise.kb = 200, indep.pairwise.ct = 30, indep.pairwise.r2 = 0.5)
#' filter.dataset(plink.files = TRUE, bed.file = "mybed.bed", bim.file = "mybim.bim", fam.file = "myfam.fam", remove.related = FALSE, geno.value = 0.5, maf.value = 0.05, mind.value = 0.1, limit.LD = TRUE, indep.pairwise.kb = 300, indep.pairwise.ct = 25, indep.pairwise.r2 = 0.4)

get_plink_path <- function(version = "plink") {
   plink <- Sys.which(version)
   if (plink == "") stop(paste(version, "not found in PATH. Ensure it's installed in /usr/local/bin"))
   return(plink)
}

run_plink <- function(args, version = "plink") {
   exec <- get_plink_path(version)
   message("Using", version, "with:", args)
   system(paste(exec, args))
}

relatedness_filter <- function(file.type, input.file, kinship.coefficient) {
   prefix <- switch(file.type,
                    VCF = paste("--vcf", input.file),
                    BCF = paste("--bcf", input.file),
                    PLINK = paste("--bfile", input.file))
   
   run_plink(paste(prefix, "--make-king triangle bin --out related"), "plink2")
   run_plink(paste(prefix, "--king-cutoff related", kinship.coefficient, "--make-bed --out unrelated"), "plink2")
   
   return("unrelated")  # name of output base from plink2
}

apply_filters <- function(basefile, geno, mind, maf) {
   run_plink(paste("--bfile", basefile,
                   "--geno", geno,
                   "--mind", mind,
                   "--maf", maf,
                   "--make-bed --out filtered"))
   return("filtered")
}

apply_ld_pruning <- function(input.base, kb, ct, r2) {
   run_plink(paste("--bfile", input.base,
                   "--indep-pairwise", kb, ct, r2,
                   "--recode vcf --out filtered.pruned"))
}

#' @export
filter.dataset <- function(input.file,
                           plink.files = FALSE,
                           bed.file = NULL,
                           bim.file = NULL,
                           fam.file = NULL,
                           remove.related = FALSE,
                           kinship.coefficient = NULL,
                           geno.value = 0.1,
                           maf.value = 0.01,
                           mind.value = 0.1,
                           limit.LD = FALSE,
                           indep.pairwise.kb = NULL,
                           indep.pairwise.ct = NULL,
                           indep.pairwise.r2 = NULL) {
   
   file.type <- detect_file_type(input.file, bed.file, bim.file, fam.file)
   
   base_input <- if (file.type == "PLINK") {
      # build --bfile argument as a single tag
      paste("--bfile", tools::file_path_sans_ext(bed.file))
   } else {
      input.file
   }
   
   if (remove.related) {
      base_input <- perform_relatedness_filter(file.type, base_input, kinship.coefficient)
   }
   
   filtered_base <- apply_filters(base_input, geno = geno.value, mind = mind.value, maf = maf.value)
   
   if (limit.LD) {
      apply_ld_pruning(filtered_base, indep.pairwise.kb, indep.pairwise.ct, indep.pairwise.r2)
   } else {
      message("LD pruning not applied.")
   }
   
   message("Dataset filtering complete.")
}