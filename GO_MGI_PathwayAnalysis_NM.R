############################################################################################################
# Fisher's test for functional enrichment in gene sets- Nick's script
############################################################################################################
options(scipen = 999)

#TERMS_PATH <- "~/Documents/Bioinformatics/GeneSets/GO_MGI/MAMMALIAN_PHENOTYPE/MOUSE/MOUSE_MAGMA_FORMAT/ANNOTATION_04_07_2018-19_32/MGI_single_gene_pheno_to_human_10_MAGMA.txt"
TERMS_PATH <- "~/PhD/pathway_analysis_annotations/ANNOTATION_05_04_2019-02_13/GENE_SETS/GO/GO_MAGMA_FORMAT/GO_ALL_PC_10-2000_MAGMA.txt"

GENE_SET_PATH <- "~/PhD/LTP_seq/pathway_analysis/prep_for_PA/output/trap_30.txt"

BACKGROUND_PATH <- "~/PhD/LTP_seq/pathway_analysis/prep_for_PA/output/background.txt"

OUT_PATH <- "~/PhD/LTP_seq/pathway_analysis/output/trap_30_pathways.txt"

# Column containing entrez gene IDs in gene set and background files
ENTREZ_COL <- 1

# Optional pathway refinement
REFINEMENT <- TRUE

# P-value significance threshold
ALPHA <- 0.05 

###############################################################################################################
# Load required packages
###############################################################################################################

if(!require(dplyr)) {
  install.packages("dplyr")
  require(dplyr)
}

if(!require(pbapply)) {
  install.packages("pbapply")
  require(pbapply)
}
library(pbapply)
###############################################################################################################
# Functions
###############################################################################################################

# Fisher test from contingency table of gene membership
fisher_test_genes <- function(term, GS, BKG) {
  
  # Create vector of genes belonging to term
  term_genes = unlist(strsplit(terms$genes[terms$term_name == term[["term_name"]]], split = " ", fixed = TRUE)) 
  
  # Remove terms already selected during refinement, if any
  if (nrow(terms_refined) > 0) {
    term_genes = remove_refined(term_genes)
  }
  
  # Genes in both the gene set and the term
  GS_term = length(intersect(term_genes, GS)) 
  
  # Genes in the gene set but not in the term
  GS_no_term = length(setdiff(GS, term_genes)) 
  
  # Genes in the term but not the gene set (all genes must exist in the background)
  term_no_GS = length(intersect(setdiff(term_genes,GS), BKG)) 
  
  # Genes in neither the gene set nor the term
  BKG_no_term_GS = length(setdiff(BKG, union(term_genes, GS))) 
  
  # Contingency table
  cont_tab = matrix(c(GS_term, GS_no_term, term_no_GS, BKG_no_term_GS), ncol = 2) 
  
  # Fisher's exact test
  fisher_out = fisher.test(cont_tab, alternative = "greater") # Fisher's exact test
  
  # Extract term, size, overlap, odds ratio and p-value
  c(term[["term_name"]], length(term_genes), length(intersect(GS, term_genes)), fisher_out$estimate, fisher_out$p.value)
}


# Remove terms already selected during refinement
remove_refined <- function(term_genes) {
  
  # Create vector of unique genes from terms already selected during refinement
  refined_genes = get_refined_genes(terms_refined)
  
  # Remove genes from current term
  setdiff(term_genes, refined_genes)
}


# Create vector of genes from data frame of terms
get_refined_genes <- function(terms_refined) {
  
  # For each term already selected for output
  unlist(apply(terms_refined, 1, function(term) {
    
    # Find genes from source file and separate
    strsplit(terms$genes[terms$term == term[["term_name"]]], split = " ")
  }))
  
}

###############################################################################################################
# Main
###############################################################################################################

# Pathway Analysis

GS <- read.table(GENE_SET_PATH, sep = "\t", header = T, stringsAsFactors = F, colClasses = "character") %>% 
  dplyr::select(ENTREZ_COL) %>% 
  unlist

BKG <- read.table(BACKGROUND_PATH, sep = "\t", header = T, stringsAsFactors = F, colClasses = "character") %>% 
  dplyr::select(ENTREZ_COL) %>% 
  unlist

terms <- read.table(TERMS_PATH, stringsAsFactors = FALSE, colClasses = "character", sep = "\t", quote = NULL, col.names = c("term_name", "genes"))

# Prep empty refinement output variable
terms_refined <- data.frame("term_name" = character(),
                            "n_term" = numeric(),
                            "n_overlap" = numeric(),
                            "odds_ratio" = numeric(),
                            "p" = numeric(), 
                            "bonferroni" = numeric())

# Compare overlap of gene set with each term using Fisher's exact test
print("Performing initial pathway analysis")
terms_enrich <- pbapply(terms, 1, fisher_test_genes, GS = GS, BKG = BKG)

# Convert to dataframe with appropriate column classes
terms_enrich <- data.frame("term_name" = as.character(terms_enrich[1,]), 
                           "n_term" = as.numeric(terms_enrich[2,]),
                           "n_overlap" = as.numeric(terms_enrich[3,]), 
                           "odds_ratio" = as.numeric(terms_enrich[4,]), 
                           "p" = as.numeric(terms_enrich[5,]),
                           stringsAsFactors = F)

# Append corrected p-values and sort
terms_enrich$bonferroni <- p.adjust(terms_enrich$p, method = "bonferroni")
terms_enrich <- arrange(terms_enrich, p)

write.table(terms_enrich, OUT_PATH, sep = "\t", row.names = FALSE, quote = FALSE)

# Refinement

if(REFINEMENT == TRUE) {
  
  # Filter out non-significant terms
  terms_sig <- filter(terms_enrich, bonferroni < ALPHA) 
  print(paste0("Finding principle terms from ", as.character(nrow(terms_sig)), " terms"))
  
  # While there are terms still significant terms needing refinement
  while (nrow(terms_sig) > 0) {
    
    # Sort terms that are still significant by odds_ratio
    terms_sig <- arrange(terms_sig, desc(odds_ratio))
    
    # Append term with highest odds ratio to output
    terms_refined <- rbind(terms_refined, terms_enrich[terms_enrich$term_name == terms_sig$term_name[1], ])
    
    # For each significant term still needing refinement, repeat Fisher's exact test
    temp_terms_enrich <- apply(terms_sig, 1, fisher_test_genes, GS = GS, BKG = BKG) 
    
    # Convert to dataframe with appropriate column classes
    temp_terms_enrich <- data.frame("term_name" = as.character(temp_terms_enrich[1,]), 
                                    "n_term" = as.numeric(temp_terms_enrich[2,]),
                                    "n_overlap" = as.numeric(temp_terms_enrich[3,]), 
                                    "odds_ratio" = as.numeric(temp_terms_enrich[4,]), 
                                    "p" = as.numeric(temp_terms_enrich[5,]),
                                    stringsAsFactors = F)
    
    # Filter to remove terms no longer significant
    terms_sig <- filter(temp_terms_enrich, p < ALPHA)
    
    print(paste0("Number of refined genes = ", length(get_refined_genes(terms_refined))))
    
  }
  
  write.table(terms_refined, paste0(dirname(OUT_PATH), "/refined_", basename(OUT_PATH)), sep = "\t", row.names = FALSE, quote = FALSE)
}


