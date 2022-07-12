#!/usr/bin/env Rscript

# Maintained by Mike Lee (Mike.Lee@nasa.gov)
# GeneLab script for generating organism ENSEMBLE annotation tables
# Example usage: Rscript build-genome-annots-tab.R MOUSE

############################################
############ pre-flight checks #############
############################################

args <- commandArgs(trailingOnly = TRUE)

currently_accepted_orgs <- c("MOUSE", "HUMAN", "ARABIDOPSIS", "FLY")

# making sure a positional argument was provided
if ( length(args) < 1 ) {
    cat("\n  One positional argument is required that specifies the target organism. Currently available include:\n")

    for ( item in currently_accepted_orgs ) {

        cat(paste0("\n        ", item))
    }

    cat("\n\n")

    quit()

} else {

    suppressWarnings(target_organism <- toupper(args[1]))

}

# making sure it is one we are prepared to handle
if ( ! target_organism %in% currently_accepted_orgs ) {

    cat(paste0("\n  '", args[1], "' isn't a valid entry.\n"))

    cat("\n  The currently available organisms include:\n")

    for ( item in currently_accepted_orgs ) {

        cat(paste0("\n        ", item))
    }

    cat("\n\n")

    quit()

}


############################################
######### setting some things up ###########
############################################

library(tidyverse)
library(STRINGdb)
library(PANTHER.db)
library(rtracklayer)


# setting primary keytype, right now, TAIR if arabidopsis, ENSEMBL if anything else
if ( target_organism == "ARABIDOPSIS" ) {

    primary_keytype <- "TAIR"

} else {

    primary_keytype <- "ENSEMBL"

}

wanted_keys_vec <- c("SYMBOL", "GENENAME", "REFSEQ", "ENTREZID")
org_tab_link <-
    "https://raw.githubusercontent.com/asaravia-butler/GeneLab_Data_Processing/master/RNAseq/organisms.csv"
ref_tab_link <-
    "https://raw.githubusercontent.com/asaravia-butler/GeneLab_Data_Processing/master/RNAseq/GeneLab_Reference_and_Annotation_Files/GL-DPPD-7101-E_ensembl_refs.csv"


############################################
######## getting down to business ##########
############################################

# setting this so any downloads won't timeout
options(timeout = 600)

# getting org taxid and org db name
organism_table <- read.csv(org_tab_link)

ref_table <- read.csv(ref_tab_link)

target_taxid <- organism_table %>%
    filter(name == target_organism) %>%
    pull(taxon)

target_org_db <- organism_table %>%
    filter(name == target_organism) %>%
    pull(annotations)

target_species_designation <- organism_table %>%
    filter(name == target_organism) %>%
    pull(species)

gtf_link <- ref_table %>%
    filter(Organism == target_species_designation) %>%
    pull(Annotation.File)

# making output prefix
base_gtf_filename <- basename(gtf_link)
base_output_name <- str_replace(base_gtf_filename, ".gtf.gz", "")

out_table_filename <- paste0(base_output_name, "-GL-annotations.tsv")
out_log_filename <- paste0(base_output_name, "-GL-build-info.txt")

# making sure output file doesn't exist already, if it does, we're exiting without overwriting
if ( file.exists(out_table_filename) ) {

    cat("\n-------------------------------------------------------------------------------------------------\n")
    cat(paste0("\n  The file that would be created, '", out_table_filename, "', exists already.\n"))
    cat(paste0("\n  We don't want to overwrite it accidentally. Move it and run this again if wanting to proceed.\n"))
    cat("\n-------------------------------------------------------------------------------------------------\n")

    quit()

}

## getting all unique gene IDs
# changing link to be https if it's currently http
gtf_link <- str_replace(gtf_link, "http:", "https:")

gtf_obj <- import(gtf_link)
unique_IDs <- gtf_obj$gene_id %>% unique()

# removing gtf object now that we don't need it cause it can soak up a lot of RAM
rm(gtf_obj)

ann.dbi <- target_org_db

# installing org db if needed
if ( ! require(ann.dbi, character.only = TRUE)) {

    BiocManager::install(ann.dbi, ask = FALSE)

}

library(ann.dbi, character.only = TRUE)

annot <- data.frame(unique_IDs)
colnames(annot) <- primary_keytype

for ( key in wanted_keys_vec ) {

    if ( key %in% columns(eval(parse(text = ann.dbi), env = .GlobalEnv))) {

        new_list <- mapIds(eval(parse(text = ann.dbi), env = .GlobalEnv), keys = unique_IDs, keytype = primary_keytype, column = key, multiVals = "list")

        # they come as lists when we accept the multiple hits, so converting to character strings here
        annot[[key]] <- sapply(new_list, paste, collapse = "|")

    }
}

# adding STRINGdb annots
string_db <- STRINGdb$new(version = "11", species = target_taxid, score_threshold = 0)
string_map <- string_db$map(annot, primary_keytype, removeUnmappedRows = FALSE, takeFirst = FALSE)

# adding a return because the stringdb stdout doesn't have one and it bothers me
cat("\n\n")

# combining if there are any with multiple string IDs
tab_with_multiple_STRINGids_combined <-
    data.frame(row.names = annot[[primary_keytype]])

for ( curr_gene_ID in row.names(tab_with_multiple_STRINGids_combined) ) {

    # getting current gene's STRING_ids and turning into
    # combined a vector
    curr_STRING_ids <- string_map %>%
        filter(!!rlang::sym(primary_keytype) == curr_gene_ID) %>%
        pull(STRING_id) %>% paste(collapse = "|")

    # adding to table
    tab_with_multiple_STRINGids_combined[curr_gene_ID, "STRING_id"] <- curr_STRING_ids

}

# moving primary_keytype back to being a column
tab_with_multiple_STRINGids_combined <-
    tab_with_multiple_STRINGids_combined %>%
    rownames_to_column(primary_keytype)

# combining string column
annot <- dplyr::left_join(annot,
                          tab_with_multiple_STRINGids_combined,
                          by = primary_keytype)

# adding GO slim annotations
pthOrganisms(PANTHER.db) <- target_organism

# since we are using ENTREZIDs to pull from the PANTHER db, and there can be multiple ENTREZIDs for a gene, we need to split them first,
# so doing this as a loop building the new GOSLIM annotation column

for ( curr_row in 1:dim(annot)[1] ) {

    curr_entry <- annot[curr_row, "ENTREZID"]

    # dealing with NAs
    if ( curr_entry == "NA" ) {

        annot[curr_row, "GOSLIM_IDS"] <- "NA"

    } else if ( ! grepl("|", curr_entry, fixed = TRUE) ) {

        # this handles if there is only one ENTREZID
        curr_GO_IDs <- mapIds(PANTHER.db, keys = curr_entry, keytype = "ENTREZ", column = "GOSLIM_ID", multiVals = "list") %>% unlist() %>% as.vector()

        # handles if none were found
        if ( is.null(curr_GO_IDs) ) {

            curr_GO_IDs <- "NA"
        }

        annot[curr_row, "GOSLIM_IDS"] <- paste(curr_GO_IDs, collapse = "|")

    } else {

        # this block handles if there are multiple ENTREZIDs, if this is the case,
        # we want to split them, get the GO IDs for all, combine them, and remove any duplicates

        # splitting them
        curr_entry_vec <- strsplit(curr_entry, "|", fixed = TRUE)

        # starting vector of current GO IDs
        curr_GO_IDs <- vector()

        # looping through, getting GO IDs, and combining them
        for ( curr_entry in curr_entry_vec ) {

            new_GO_IDs <- mapIds(PANTHER.db, keys = curr_entry, keytype = "ENTREZ", column = "GOSLIM_ID", multiVals = "list") %>% unlist() %>% as.vector()

            # adding to building vector of all GO IDs
            curr_GO_IDs <- c(curr_GO_IDs, new_GO_IDs)

        }

        # removing any dups
        curr_GO_IDs <- unique(curr_GO_IDs)

        # handles if none were found
        if ( length(curr_GO_IDs) == 0 ) {

            curr_GO_IDs <- "NA"
        }

        # adding to table
        annot[curr_row, "GOSLIM_IDS"] <- paste(curr_GO_IDs, collapse = "|")

    }

}


############################################
##### writing out annots and run info ######
############################################

# ordering first
annot <- annot %>% arrange(.[[1]])
# writing out annotations table
write.table(annot, out_table_filename, sep = "\t", quote = FALSE, row.names = FALSE)


# getting date of run
date_generated <- format(Sys.time(), "%d-%B-%Y")

# writing out building info
writeLines(paste(c("Build done on:\n    ", date_generated), collapse = ""), out_log_filename)
write(paste(c("\nUsed gtf file:\n    ", gtf_link), collapse = ""), out_log_filename, append = TRUE)
write(paste(c("\nUsed ", ann.dbi, " version:\n    ", packageVersion(ann.dbi) %>% as.character()), collapse = ""), out_log_filename, append = TRUE)
write(paste(c("\nUsed STRINGdb version:\n    ", packageVersion("STRINGdb") %>% as.character()), collapse = ""), out_log_filename, append = TRUE)
write(paste(c("\nUsed PANTHER.db version:\n    ", packageVersion("PANTHER.db") %>% as.character()), collapse = ""), out_log_filename, append = TRUE)

write("\n\nAll session info:\n", out_log_filename, append = TRUE)
write(capture.output(sessionInfo()), out_log_filename, append = TRUE)
