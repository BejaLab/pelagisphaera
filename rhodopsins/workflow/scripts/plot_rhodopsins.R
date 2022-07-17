
library(dplyr)
library(treeio)
library(ggtree)
library(phangorn)
library(bioformatr)
library(seqinr)
library(tidyr)
library(taxize)
library(tidyverse)
library(castor)
library(ggnewscale)

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots = list(input = "list", output = "list", params = "list"))
    snakemake <- Snakemake(
            input = list(
            treefile = "analysis/iqtree/combined.treefile",
            metadata = "metadata/metadata.txt",
            uniparc_tab = "input/P_uniparc.tab",
            uniprot_tab = "input/P_uniprot.tab",
            mafft = "analysis/fasta/combined_100.cdhit.mafft",
            taxonomy = "metadata/gtdbtk_taxonomy.tsv",
            clstr60  = "analysis/fasta/combined_60.cdhit.clstr",
            clstr70  = "analysis/fasta/combined_70.cdhit.clstr",
            clstr80  = "analysis/fasta/combined_80.cdhit.clstr",
            clstr90  = "analysis/fasta/combined_90.cdhit.clstr",
            clstr100 = "analysis/fasta/combined_100.cdhit.clstr",
            taxize   = "analysis/taxize/taxize.tsv",
            outgroups = "input/outgroups.fasta"
        ),
        output = list(
            "output/rhodopsins.svg"
        )
    )
}
with(snakemake@input, {
    uniparc_tab    <<- uniparc_tab
    uniprot_tab    <<- uniprot_tab
    metadata_file  <<- metadata
    mafft_file     <<- mafft
    taxonomy_file  <<- taxonomy
    clstr60_file   <<- clstr60
    clstr70_file   <<- clstr70
    clstr80_file   <<- clstr80
    clstr90_file   <<- clstr90
    clstr100_file  <<- clstr100
    treefile       <<- treefile
    taxize_file    <<- taxize
    outgroup_file  <<- outgroups
})
output_file <- unlist(snakemake@output)

read.cdhit.clstr <- function(fname) {
    data.fields <- c("E.Value", "Aln", "Identity")
    read.table(fname, sep = "\t", comment.char = "", quote = "", fill = T, stringsAsFactors = F, col.names = c("Col1", "Col2")) %>%
        separate(Col1, into = c("Seq.Num", "Cluster"), sep = " ", fill = "right") %>%
        fill(Cluster) %>%
        filter(!grepl(">", Seq.Num)) %>%
        separate(Col2, into = c("Seq.Len", "Col2"), sep = "aa, >") %>%
        extract(Col2, into = c("Seq.Name", "Is.Representative", "Col2"), regex = "(.*?)[.]{3} ([*]|at) ?(.*)") %>%
        mutate(Is.Representative = Is.Representative == "*", Col2 = ifelse(Is.Representative, "100%", Col2)) %>%
        group_by(Cluster) %>%
        mutate(Representative = Seq.Name[which(Is.Representative)]) %>%
        separate_rows(Col2, sep = ",") %>%
        separate(Col2, into = data.fields, sep = "/", fill = "left", convert = T) %>%
        mutate(Identity = sub("%", "", Identity) %>% as.numeric) %>%
        group_by(Seq.Name) %>%
        mutate(level.rank = paste0(".", 1:n() - 1), level.rank = ifelse(level.rank == ".0", "", level.rank)) %>%
        pivot_wider(names_from = level.rank, values_from = data.fields, names_sep = "") %>%
        ungroup
}

metadata <- read.table(metadata_file, sep = "\t", header = T)

uniparc <- read.table(uniparc_tab, sep = "\t", header = T) %>%
    select(Entry, Organism = Organisms, Organism.ID = Organism.IDs) %>%
    separate_rows(Organism, Organism.ID, sep = "; ", convert = T)
uniprot <- read.table(uniprot_tab, sep = "\t", header = T) %>%
    select(Entry, Organism, Organism.ID)
unitab <- bind_rows(uniparc, uniprot)

ncbi.tax <- read.table(taxize_file, sep = "\t", col.names = c("Organism.ID", "taxgroup"))

metadata <- read.table(metadata_file, sep = "\t", header = T, na.string = "") %>%
    arrange(is.na(Known)) %>%
    distinct(Entry, .keep_all = T)

#rhodopsins <- read.table("Verrucomicrobia_all.txt", header = T) %>%
#    distinct(label = Representative, family, subfamily)

motif <- c(97, 101, 108) %>% as.character
GPR <- read.fasta(mafft_file, seqtype = "AA") %>%
    lapply(c) %>%
    lapply(as.data.frame) %>%
    lapply(rowid_to_column, "pos") %>%
    bind_rows(.id = "Entry") %>%
    spread(pos, 3) %>%
    select(Entry, D85 = !!motif[1], T89 = !!motif[2], D96 = !!motif[3])

verruco <- read.table(taxonomy_file, col.names = c("assembly", "gtdb"), sep = "\t") %>%
    filter(grepl("p__Verrucomicrobiota", gtdb)) %>%
    pull(assembly)

clstr60  <- read.cdhit.clstr(clstr60_file) %>%
    select(Representative_70 = Seq.Name, Representative_60 = Representative)
clstr70  <- read.cdhit.clstr(clstr70_file) %>%
    select(Representative_80 = Seq.Name, Representative_70 = Representative)
clstr80  <- read.cdhit.clstr(clstr80_file) %>%
    select(Representative_90 = Seq.Name, Representative_80 = Representative)
clstr90  <- read.cdhit.clstr(clstr90_file) %>%
    select(Representative_100 = Seq.Name, Representative_90 = Representative)
clstr100 <- read.cdhit.clstr(clstr100_file) %>%
    select(Entry = Seq.Name, Representative_100 = Representative) %>%
    left_join(clstr90, by = "Representative_100") %>%
    left_join(clstr80, by = "Representative_90") %>%
    left_join(clstr70, by = "Representative_80") %>%
    left_join(clstr60, by = "Representative_70") %>%
    left_join(unitab, by = "Entry") %>%
    left_join(ncbi.tax, by = "Organism.ID") %>%
    left_join(GPR, by = "Entry") %>%
    left_join(metadata, by = "Entry") %>%
    filter(!is.na(Representative_60)) %>%
    separate(Entry, into = c("assembly", "locus_tag"), sep = "@", fill = "right", remove = F) %>%
    mutate(taxgroup = ifelse(!is.na(assembly) & assembly %in% verruco, "Verrucomicrobia", ifelse(taxgroup == "Verrucomicrobia", NA, taxgroup))) %>%
    group_by(Representative_60, taxgroup) %>%
    mutate(taxgroup_n = sum(!is.na(taxgroup))) %>%
    ungroup %>%
    arrange(-taxgroup_n) %>%
    group_by(Representative_60) %>%
    mutate(Known = first(na.omit(Known))) %>%
    mutate(Symbol = paste(unique(na.omit(Symbol)), collapse = ", ")) %>%
    mutate(Family = first(na.omit(Family))) %>%
    ungroup
clstr <- group_by(clstr100, Representative_60, D85, T89, D96) %>%
    mutate(n_motif = n()) %>%
    ungroup %>%
    arrange(is.na(D85), -n_motif) %>%
          distinct(label = Representative_60, .keep_all = T) %>%
    # left_join(rhodopsins,  by = "label") %>%
    mutate(label = sub("@", "_", label))

outgroups <- names(read.fasta(outgroup_file))

tree <- read.tree(treefile) %>%
    # midpoint %>%
    ape::root(outgroups, edgelabel = T, resolve.root = T) %>%
    drop.tip(outgroups) %>%
    as_tibble %>%
    mutate(support = suppressWarnings(as.numeric(label))) %>%
    left_join(clstr, by = "label") %>%
    group_by(taxgroup) %>%
    mutate(taxgroup_n = n()) %>%
    ungroup %>%
    mutate(taxgroup = ifelse(taxgroup_n < 3, NA, taxgroup)) %>%
    `class<-`(c("tbl_tree", "data.frame"))

tree.phylo <- as.treedata(tree)@phylo

families <- filter(tree, ! node %in% parent) %>%
    pull(Family) %>%
    as.factor
hsp <- hsp_max_parsimony(tree.phylo, as.numeric(families)) %>%
        `$`("likelihoods") %>%
    data.frame(check.names = F) %>%
    mutate(node = row_number()) %>%
        gather(hsp, prob, -node) %>%
        filter(prob > 0.9) %>%
        mutate(hsp = levels(families)[as.numeric(hsp)])
tree <- left_join(tree, hsp, by = "node") %>%
    mutate(hsp.parent = .[match(parent, node),"hsp"]) %>%
    mutate(hsp = ifelse(!is.na(hsp.parent) & hsp == hsp.parent, NA, hsp))

clustalx <- c(
    A = "BLUE",
    I = "BLUE",
    L = "BLUE",
    M = "BLUE",
    F = "BLUE",
    W = "BLUE",
    V = "BLUE",
    C = "BLUE",
    K = "RED",
    R = "RED",
    E = "MAGENTA",
    D = "MAGENTA",
    N = "GREEN",
    Q = "GREEN",
    S = "GREEN",
    T = "GREEN",
    G = "ORANGE",
    P = "YELLOW",
    H = "CYAN",
    Y = "CYAN"
)

# colors.families <- c(XR_E3 = "#44aa00", PR_FW = "#00aad4", PR_Methylacidiphilales = "#008080ff", PR_F3 = "#0000ff", P5 = "#800080")
colors.pumps    <- c(proton = "red", sodium = "magenta", chloride = "yellow4")
verruco.clades <- c(
        singleton                   = "darkgray",
        `other rhodopsins`          = "gray",
        `Clade_P5-A0A350AVH5`       = "#800080",
        `Clade_P5-PWZW01000279_GM_1972` = "#800080",
        `PR-DLQG01000084_GM_551`    = "#0000ff",
        `PR-VFJU01000426_GM_2367`   = "#00aad4",
        `PR-CAJAHL010000291_GM_2019`= "#8b93ff",
        `XR-F7C95_05130`            = "#44aa00",
        `XR-CAIZMQ010000037_GM_448` = "#00ff00",
        `Clade_P4-A0A2D5AQR4` = "#ff00ff",
        `Clade_P4-A0A2E5DYM8` = "#ff00ff",
        `Clade_P4-A0A2E8M222` = "#ff00ff",
        `Clade_P4-CAIVXP010000160_GM_1042` = "#ff00ff"
)

p <- ggtree(as.treedata(tree), aes(color = taxgroup), layout = "circular") +
    #geom_tiplab(aes(subset = taxgroup == "Verrucomicrobia", label = label)) +
    geom_treescale() +
    geom_tiplab(mapping = aes(subset = !is.na(Known), label = Known), color = "red", size = 2) +
    geom_point2(aes(subset = !is.na(support) & support >= 95), shape = 15, color = "darkgray") +
    new_scale_color() +
    # geom_point2(aes(subset = !is.na(subfamily), color = subfamily)) +
    scale_color_manual(values = verruco.clades) +
    # geom_tiplab2(aes(label = label)) +
    geom_tiplab(aes(subset = !is.na(Symbol), label = Symbol), offset = 0.1) +
    new_scale_color() +
    geom_text2(aes(label = D85, color = D85, angle = angle - 90, x = 5.8), size = 2) +
    geom_text2(aes(label = T89, color = T89, angle = angle - 90, x = 6.0), size = 2) +
    geom_text2(aes(label = D96, color = D96, angle = angle - 90, x = 6.2), size = 2) +
    scale_color_manual(values = clustalx) +
    geom_highlight(mapping = aes(subset = !is.na(hsp), fill = hsp), alpha = 0.1)
ggsave(output_file, p, width = 10, height = 10)
