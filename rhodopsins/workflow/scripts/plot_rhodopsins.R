
library(dplyr)
library(treeio)
library(ggtree)
library(phangorn)
library(seqinr)
library(tidyr)
library(taxize)
library(tidyverse)
library(castor)
library(ggnewscale)

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
    family_colors_file <<- family_colors
    rhodopsins_file <<- rhodopsins
})
output_file <- unlist(snakemake@output)

rhodopsins <- read.table(rhodopsins_file, header = T) %>%
    mutate(Entry = paste(label, locus_tag, sep = "@"))

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

uniparc <- read.table(uniparc_tab, sep = "\t", header = T) %>%
    select(Entry, Organism = Organisms, Organism.ID = Organism.IDs) %>%
    separate_rows(Organism, Organism.ID, sep = "; ", convert = T)
uniprot <- read.table(uniprot_tab, sep = "\t", header = T) %>%
    select(Entry, Organism, Organism.ID)
unitab <- bind_rows(uniparc, uniprot)

ncbi.tax <- read.table(taxize_file, sep = "\t", col.names = c("Organism.ID", "taxgroup"))

metadata <- read.table(metadata_file, sep = "\t", header = T, na.string = "", fill = T) %>%
    arrange(is.na(Activity)) %>%
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

verruco <- read.table(taxonomy_file, col.names = c("label1", "label2", "gtdb"), sep = "\t") %>%
    filter(grepl("p__Verrucomicrobiota", gtdb)) %>%
    select(label1, label2) %>%
    as.matrix %>% c %>%
    unique
colors.families <- read.table(family_colors_file, sep = "\t", comment.char = "") %>%
    with(setNames(V2, V1))

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
    left_join(rhodopsins, by = "Entry") %>%
    filter(!is.na(Representative_60)) %>%
    separate(Entry, into = c("assembly", "locus_tag"), sep = "@", fill = "right", remove = F) %>%
    mutate(taxgroup = ifelse(!is.na(assembly) & assembly %in% verruco, "Verrucomicrobia", ifelse(taxgroup == "Verrucomicrobia", NA, taxgroup))) %>%
    group_by(Representative_60, taxgroup) %>%
    mutate(taxgroup_n = sum(!is.na(taxgroup))) %>%
    ungroup %>%
    arrange(-taxgroup_n) %>%
    group_by(Representative_60) %>%
    mutate(Activity = first(na.omit(Activity))) %>%
    mutate(Symbol = paste(unique(na.omit(Symbol)), collapse = ", ")) %>%
    mutate(Family = first(na.omit(Family))) %>%
    ungroup
clstr <- group_by(clstr100, Representative_60, D85, T89, D96) %>%
    mutate(n_motif = n(), subfamily = first(na.omit(subfamily))) %>%
    ungroup %>%
    arrange(is.na(D85), -n_motif) %>%
    distinct(label = Representative_60, .keep_all = T) %>%
    mutate(label = sub("@", "_", label))

missing_fams <- filter(clstr, !is.na(subfamily)) %>%
    filter(! subfamily %in% names(colors.families)) %>%
    distinct(subfamily)
if (nrow(missing_fams) > 0) {
    write(paste("The following families have unassigned colors:", paste(pull(missing_fams), collapse = ", ")), stderr())
    q()
}
missing_fams <- data.frame(fam = names(colors.families)) %>%
    filter(! fam %in% clstr$subfamily, !grepl(" ", fam)) %>%
    distinct(fam)
if (nrow(missing_fams) > 0) {
    write(paste("The following families have were not found in the data:", paste(pull(missing_fams), collapse = ", ")), stderr())
    q()
}

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
    mutate(hsp = ifelse(!is.na(hsp.parent) & hsp == hsp.parent, NA, hsp)) %>%
    mutate(Ion = case_when(Activity == "chloride pump" ~ "Cl⁻", Activity == "sodium pump" ~ "Na⁺", Activity == "proton pump" ~ "H⁺", T ~ NA_character_))

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
p <- ggtree(as.treedata(tree), aes(color = taxgroup), layout = "circular") +
    #geom_tiplab(aes(subset = taxgroup == "Verrucomicrobia", label = label)) +
    geom_treescale() +
    geom_tiplab(mapping = aes(subset = !is.na(Ion), label = Ion), color = "red", size = 2) +
    geom_point2(aes(subset = !is.na(support) & support >= 95), shape = 15, color = "darkgray") +
    new_scale_color() +
    geom_point2(aes(subset = !is.na(subfamily), color = subfamily)) +
    scale_color_manual(values = colors.families) +
    # geom_tiplab2(aes(label = label)) +
    geom_tiplab(aes(subset = !is.na(Symbol) & !is.na(Ion), label = Symbol), offset = 0.3) +
    new_scale_color() +
    geom_text2(aes(label = D85, color = D85, angle = angle - 90, x = 5.8), size = 2) +
    geom_text2(aes(label = T89, color = T89, angle = angle - 90, x = 6.0), size = 2) +
    geom_text2(aes(label = D96, color = D96, angle = angle - 90, x = 6.2), size = 2) +
    scale_color_manual(values = clustalx) +
    geom_highlight(mapping = aes(subset = !is.na(hsp), fill = hsp), alpha = 0.1)
ggsave(output_file, p, width = 10, height = 10)
