library(dplyr)
library(tidyr)
library(tibble)
library(seqinr)

with(snakemake@input, {
    usearch_file <<- usearch
    blastp_file  <<- blastp
})
output_file <- unlist(snakemake@output)

read.outfmt6 <- function(fname, col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) {
    read.table(fname, sep = "\t", col.names = col.names)
}

blast <- read.outfmt6(usearch_file) %>%
    select(Representative = sseqid, Seq.Name = qseqid)
read.outfmt6(blastp_file, col.names = c("Seq.Name", "pident", "stitle")) %>%
    extract(stitle, into = c("superfamily", "family"), regex = "/([\\w?]+):(\\w+)/") %>%
    left_join(blast, by = "Seq.Name") %>%
    separate(Seq.Name, into = c("label", "locus_tag"), sep = "@") %>%
    separate(Representative, into = c("Representative_label", "Representative_locus_tag"), sep = "@", fill = "left", remove = F) %>%
    group_by(family, Representative) %>%
    mutate(n_fam_repr = n()) %>%
    group_by(family) %>%
    mutate(n_fam = n()) %>%
    mutate(subfamily = case_when(is.na(Representative) ~ "other rhodopsins", n_fam < 3 | n_fam_repr == 1 ~ "singleton", T ~ paste(family, Representative_locus_tag, sep = "-"))) %>%
    ungroup %>%
    arrange(-n_fam_repr) %>%
    write.table(file = output_file, row.names = F)
