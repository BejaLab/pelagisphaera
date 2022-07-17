
library(dplyr)
library(tidyr)
library(taxize)

input_files <- unlist(snakemake@input)
output_file <- unlist(snakemake@output)

lapply(input_files, read.table, sep = "\t", header = T) %>%
    lapply(rename, any_of(c(Organism.ID = "Organism.IDs"))) %>%
    lapply(separate_rows, "Organism.ID", sep = "; ", convert = T) %>%
    bind_rows %>%
    distinct(Organism.ID) %>%
    pull %>%
    classification(db = "ncbi") %>%
    `[`(!is.na(.)) %>%
    bind_rows(.id = "Organism.ID") %>%
    mutate(Organism.ID = as.numeric(Organism.ID)) %>%
    select(Organism.ID, name, rank) %>%
    distinct(Organism.ID, rank, .keep_all = T) %>%
    spread(rank, name) %>%
    mutate(taxgroup = ifelse(superkingdom == "Bacteria", phylum, superkingdom)) %>%
    select(Organism.ID, taxgroup) %>%
    write.table(file = output_file, sep = "\t", row.names = F, col.names = F)
