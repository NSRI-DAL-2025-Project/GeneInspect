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
#' @export

filter.dataset <- function(
      input.file = input, 
      plink.files = FALSE,
      bed.file = bed,
      bim.file = bim,
      fam.file = fam,
      remove.related = FALSE,
      kinship.coefficient = NULL,
      geno.value = geno,
      maf.value = maf,
      mind.value = mind,
      limit.LD = FALSE,
      indep.pairwise.kb = NULL,
      indep.pairwise.ct = NULL,
      indep.pairwise.r2 = NULL){
   
   plink_exec <- function(PLINKoptions = "") system(paste("src/plink.exe",PLINKoptions)) #specify path to plink1.9
   plink_exec2 <- function(PLINKoptions = "") system(paste("src/plink2.exe",PLINKoptions)) #specify path to plink2.0
   
   
   if(plink.files == TRUE){
      input <- sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file) # read in plink files
      
      if(remove.related == TRUE){
         # checking related samples
         plink_exec2(stringr::str_c("--bfile ", 
                                    input, 
                                    " --make-king triangle bin --out related"))
         
         plink_exec2(stringr::str_c("--bfile ",
                                    input, 
                                    " --king-cutoff related ", 
                                    kinship.coefficient, 
                                    " --make-bed --out unrelated"))
         
         plink_exec(stringr::str_c(
            "--bfile unrelated --geno ",
            geno.value,
            " --mind ", 
            mind.value,
            " --maf ",
            maf.value,
            " --make-bed --out filtered"))
         
         if(limit.LD == FALSE){
            print("No LD pruning.")
         } else if (limit.LD == TRUE){
            plink_exec(stringr::str_c(
               "--bfile filtered --indep-pairwise ",
               indep.pairwise.kb,
               indep.pairwise.ct, 
               indep.pairwise.r2,
               " --recode vcf --out filtered.pruned"))
         }
         
      } else if(remove.related == FALSE){
         
         plink_exec(stringr::str_c(
            "--bfile ",
            input,
            " --geno ",
            geno.value,
            " --mind ", 
            mind.value,
            " --maf ",
            maf.value,
            " --make-bed --out filtered"))
         
         if(limit.LD == FALSE){
            print("No LD pruning.")
         } else if (limit.LD == TRUE){
            plink_exec(stringr::str_c(
               "--bfile filtered --indep-pairwise ",
               indep.pairwise.kb,
               indep.pairwise.ct, 
               indep.pairwise.r2,
               " --recode vcf --out filtered.pruned"))
         }
         
      }
      
      
   } else if(plink.files == FALSE){
      
      if(remove.related == TRUE){
         if(tools::file_ext(input.file) == "vcf"){
            plink_exec2(stringr::str_c(
               "--vcf ", input.file, 
               " --make-king triangle bin --out related"))
            
            plink_exec2(stringr::str_c("--vcf ", input.file, 
                                       " --king-cutoff related ", 
                                       kinship.coefficient, 
                                       " --make-bed --out unrelated"))
            
            # for bcf input files
         } else if(tools::file_ext(input.file) == "bcf"){
            plink_exec2(stringr::str_c(
               "--bcf ", input.file, 
               " --make-king triangle bin --out related"))
            
            plink_exec2(stringr::str_c("--bcf ", input.file, 
                                       " --king-cutoff related ", 
                                       kinship.coefficient, 
                                       " --make-bed --out unrelated"))
         }
         
         plink_exec(stringr::str_c(
            "--bfile unrelated --geno ",
            geno.value,
            " --mind ", 
            mind.value,
            " --maf ",
            maf.value,
            " --make-bed --out filtered"))
         
         if(limit.LD == FALSE){
            print("No LD pruning.")
         } else if (limit.LD == TRUE){
            plink_exec(stringr::str_c(
               "--bfile filtered --indep-pairwise ",
               indep.pairwise.kb,
               indep.pairwise.ct, 
               indep.pairwise.r2,
               " --recode vcf --out filtered.pruned"))
         }
         
      } else if(remove.related == FALSE){
         if(tools::file_ext(input.file) == "vcf"){
            
            plink_exec(stringr::str_c(
               "--vcf ",
               input.file,
               " --geno ",
               geno.value,
               " --mind ", 
               mind.value,
               " --maf ",
               maf.value,
               " --make-bed --out filtered"))
            
            if(limit.LD == FALSE){
               print("No LD pruning.")
            } else if (limit.LD == TRUE){
               plink_exec(stringr::str_c(
                  "--bfile filtered --indep-pairwise ",
                  indep.pairwise.kb,
                  indep.pairwise.ct, 
                  indep.pairwise.r2,
                  " --recode vcf --out filtered.pruned"))
            }
            
         } else if(tools::file_ext(input.file) == "bcf"){
            
            plink_exec(stringr::str_c(
               "--bcf ",
               input.file,
               " --geno ",
               geno.value,
               " --mind ", 
               mind.value,
               " --maf ",
               maf.value,
               " --make-bed --out filtered"))
            
            if(limit.LD == FALSE){
               print("No LD pruning.")
            } else if (limit.LD == TRUE){
               plink_exec(stringr::str_c(
                  "--bfile filtered --indep-pairwise ",
                  indep.pairwise.kb,
                  indep.pairwise.ct, 
                  indep.pairwise.r2,
                  " --recode vcf --out filtered.pruned"))
            }
         }
      }  
      
   } 
}
