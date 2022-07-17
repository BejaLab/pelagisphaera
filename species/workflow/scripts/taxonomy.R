library(tools)
library(gplots)
library(Biostrings)
library(readr)
library(stringr)
library(treeio)
library(ggtree)
library(tibble)
library(tidyr)
library(dplyr)
library(seqinr)
library(ggstance)

with(snakemake@input, {
    ssu_long_file   <<- ssu_long
    ssu_all_file    <<- ssu_all
    cluster_file    <<- clusters
    gtdbtk_file     <<- gtdbtk
    rhodopsins_file <<- rhodopsins
    assemblies_file <<- assemblies
    tree_file       <<- tree
    drep_file       <<- drep
})
with(snakemake@params, {
    outgroups <<- outgroups
})

ssu_long <- read.fasta(ssu_long_file, as.string = T, forceDNAtolower = F)
ssu_all  <- read.fasta(ssu_all_file,  as.string = T, forceDNAtolower = F)

clusters <- read.table(cluster_file, header = T) %>%
    rename(label = 1)
gtdbtk <- read.table(gtdbtk_file, sep = "\t", col.names = c("label", "taxon")) %>%
    separate_rows(taxon, sep = ";") %>%
    filter(grepl("g__", taxon))
rhodopsins <- read.table(rhodopsins_file, header = T) %>%
    distinct(label, .keep_all = T)

assemblies <- read.table(assemblies_file, sep = "\t", col.names = c("label", "organism")) %>%
    left_join(gtdbtk,     by = "label") %>%
    left_join(clusters,   by = "label") %>%
    left_join(rhodopsins, by = "label") %>%
    replace_na(list(taxon = ""))

organisms <- assemblies %>%
    mutate(new.ID = paste(label, taxon, organism, sep = "@")) %>%
    with(setNames(new.ID, label))

tree <- read.tree(tree_file) %>%
    ape::root(outgroups) %>%
    as_tibble %>%
    left_join(assemblies, by = "label") %>%
    mutate(support = ifelse(node %in% parent, label, NA)) %>%
    separate(support, into = c("SH.aLRT", "UFBoot"), sep = "/") %>%
    mutate(SH.aLRT = recode(SH.aLRT, `100` = "✱"), UFBoot = recode(UFBoot, `100` = "✱")) %>%
    `class<-`(c("tbl_tree", "tbl_df", "data.frame"))
family.colors <- c(PR = "#0000ff", XR = "#44aa00")

metadata <- filter(tree, ! node %in% parent) %>%
    select(id = label, Complete = Completeness)
p <- ggtree(as.treedata(tree)) +
    geom_tiplab(aes(label = label),    align = T) +
    geom_tiplab(aes(label = organism), align = T, offset = 0.1) +
    geom_text2(aes(subset = !is.na(SH.aLRT), label = SH.aLRT), hjust = 1, nudge_y =  0.25) +
    geom_text2(aes(subset = !is.na(UFBoot),  label = UFBoot),  hjust = 1, nudge_y = -0.25) +
    geom_tippoint(aes(subset = !is.na(family), color = family), size = 3) +
    scale_color_manual(values = family.colors) +
    geom_treescale(width = 0.1)
p <- facet_plot(p, data = metadata, panel = "Completeness", geom = geom_barh, mapping = aes(x = Complete, fill = Complete), stat = "identity")

ggsave(snakemake@output$phylogeny, p, width = 10, height = 6)

dists <- list()

EDNAFULL.base <- nucleotideSubstitutionMatrix(match = 5, mismatch = -4, baseOnly = TRUE)
PID2 <- Vectorize(function(x,y) pid(pairwiseAlignment(x, y, gapOpening = 16, gapExtension = 4, substitutionMatrix = EDNAFULL.base), "PID2"))

dists$SSU       <- outer(ssu_long, ssu_long, PID2) %>% `/`(100)
dists$SSU_short <- outer(ssu_all,  ssu_all,  PID2) %>% `/`(100)

ndb <- read.table(drep_file, header = T, sep = ",") %>%
    rename(A = 1, B = 2, ANI = ani, AF = alignment_coverage) %>%
    mutate(A = file_path_sans_ext(A), B = file_path_sans_ext(B))

dists$ANI <- xtabs(ANI ~ A + B, ndb)
dists$AF  <- xtabs(AF  ~ A + B, ndb)

cellnotes <- list(
    AF  = function(x) ifelse(x >= 0.6, "#", ifelse(x > 0.444,  "***", ifelse(x > 0.345,  "**", ifelse(x > 0.206,  "*", "")))),
    ANI = function(x) ifelse(x >= 0.965, "#", ifelse(x > 0.7656, "***", ifelse(x > 0.7311, "**", ifelse(x > 0.7085, "*", "")))),
    SSU = function(x) ifelse(x > 0.964,  "***", ifelse(x > 0.945,  "**", "")),
    SSU_short = function(x) ifelse(x > 0.945,  "***", "")
)

invisible(lapply(names(dists), function(V) {
    out_file <- snakemake@output[[V]]
    pdf(out_file)
    val <- dists[[V]]
    rownames(val) <- recode(rownames(val), !!!organisms)
    heatmap.2(val, col = bluered(100), cellnote = cellnotes[[V]](val), margins = c(10,10), trace = "none", distfun = function(x) as.dist(1 - x), hclustfun = function(x) hclust(x, method = "average"), vline = 0.7)
    dev.off()
}))
