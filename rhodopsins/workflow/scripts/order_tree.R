library(dplyr)
library(tidyr)
library(ggtree)
library(treeio)
library(tools)
library(phangorn)
library(ggplot2)
library(tibble)
library(ggnewscale)
library(castor)
library(ggrepel)
# library(phytools)

options(expressions = 10000)

with(snakemake@input, {
    biosamples_file <<- biosamples
    phylophlan_tree <<- phylophlan
    gtdbtk_tree     <<- gtdbtk_tree
    gtdbtk_taxonomy <<- gtdbtk_taxonomy
    dRep_dir        <<- dRep
    rhodopsins_file <<- rhodopsins
    family_colors_file <<- family_colors
    habitat_colors_file <<- habitat_colors
})
with(snakemake@params, {
    order_name <<- order
})
output_file <- unlist(snakemake@output)

phylophlan.phylo <- read.tree(phylophlan_tree)

gtdbtk   <- read.tree(gtdbtk_tree)
root.node <- which(grepl(order_name, gtdbtk$node.label)) + length(gtdbtk$tip.label)
outgroup <- Descendants(gtdbtk, root.node, "children") %>%
    Descendants(gtdbtk, ., "tips") %>%
    `[[`(2) %>% `[`(gtdbtk$tip.label, .) %>%
    gsub("^(GB|RS)_", "", .) %>%
    cat(sep = "\n", file = "Opitutales.")

taxonomy <- read.table(gtdbtk_taxonomy, col.names = c("label1", "label2", "classification"), header = F, sep = "\t")
taxonomy <- bind_rows(
         select(taxonomy, label = label1, classification),
         select(taxonomy, label = label2, classification)
    ) %>%
    distinct(label, classification) %>%
    separate_rows(classification, sep = ";") %>%
    separate(classification, into = c("rank","taxon"), sep = "__") %>%
    spread(rank, taxon) %>%
    select(label, o, f, g, s)

cdb <- file.path(dRep_dir, "data_tables/Cdb.csv") %>%
    read.table(sep = ",", header = T) %>%
    mutate(label = file_path_sans_ext(genome)) %>%
    mutate(on.tree = label %in% phylophlan.phylo$tip.label) %>%
    group_by(secondary_cluster) %>%
    filter(any(on.tree)) %>%
    mutate(cluster_label = label[on.tree])

metadata <- read.table(biosamples_file, header = T) %>%
    left_join(cdb, by = "label") %>%
    filter(!is.na(cluster_label)) %>%
    mutate(cultured_organism = ifelse(is.cultured, Organism, NA)) %>%
    group_by(cluster_label) %>%
    mutate(cultured_organism = first(na.omit(cultured_organism)))

rhodopsins <- read.table(rhodopsins_file, header = T) %>%
    left_join(cdb, by = "label") %>%
    filter(!is.na(cluster_label)) %>%
    group_by(cluster_label, subfamily) %>%
    mutate(n = n()) %>%
    arrange(-n) %>%
    ungroup %>%
    distinct(cluster_label, .keep_all = T) %>%
    select(label = cluster_label, family, subfamily)

fill_rank <- function(.data, rank.new, tree) {
    group_by(.data, !!rank.new) %>%
        mutate(mrca = getMRCA(tree, node[is.tip]) %>% replace(is.null(.), NA), mrca = ifelse(is.na(mrca), node, mrca)) %>%
        group_by(mrca) %>%
        mutate(is.mrca = node == mrca & sum(is.tip) > 0,
            rank = ifelse(is.mrca, rank.new, rank),
            label = ifelse(is.mrca, first(na.omit(!!rank.new)), label),
            !!rank.new := ifelse(is.mrca, label, !!rank.new))
}

phylophlan <- as_tibble(phylophlan.phylo) %>%
    mutate(is.tip = ! node %in% parent, is.root = node == parent, rank = ifelse(is.tip, "_", NA)) %>%
    left_join(rhodopsins, by = "label") %>%
    left_join(taxonomy, by = "label") %>%
    left_join(metadata, by = "label") %>%
    mutate(support = ifelse(node %in% parent, label, NA)) %>%
    separate(support, into = c("SH_aLRT", "ufboot"), sep = "/", convert = T) %>%
    # fill_rank("f", phylophlan.phylo) %>%
    mutate(fg = paste(f, g)) %>%
    mutate(f = ifelse(f == "", label, f)) %>%
    mutate(g = ifelse(g == "", label, g))

gens <- filter(phylophlan, is.tip) %>%
    pull(g) %>%
    as.factor
asr_gens <- asr_max_parsimony(phylophlan.phylo, as.numeric(gens)) %>%
    `$`("ancestral_likelihoods") %>%
    data.frame(check.names = F) %>%
    mutate(node = row_number() + Ntip(phylophlan.phylo)) %>%
    gather(g.asr, g.prob, -node) %>%
    filter(g.prob > 0.9) %>%
    mutate(g.asr = levels(gens)[as.numeric(g.asr)])

fams <- filter(phylophlan, is.tip) %>%
    pull(f) %>%
    as.factor
asr_fams <- asr_max_parsimony(phylophlan.phylo, as.numeric(fams), weight_by_scenarios = F) %>%
    `$`("ancestral_likelihoods") %>%
    data.frame(check.names = F) %>%
    mutate(node = row_number() + Ntip(phylophlan.phylo)) %>%
    gather(f.asr, f.prob, -node) %>%
    filter(f.prob > 0.9) %>%
    mutate(f.asr = levels(fams)[as.numeric(f.asr)])

phylophlan.tree <- phylophlan %>%
    mutate(branch.length = ifelse(branch.length < 0, 0, branch.length)) %>%
    left_join(asr_fams, by = "node") %>%
    left_join(asr_gens, by = "node") %>%
    mutate(f.asr.parent = .[match(parent, node),"f.asr"]) %>%
    mutate(g.asr.parent = .[match(parent, node),"g.asr"]) %>%
    mutate(f.asr = ifelse(!is.na(f.asr.parent) & f.asr == f.asr.parent, NA, f.asr)) %>%
    mutate(g.asr = ifelse(!is.na(g.asr.parent) & g.asr == g.asr.parent, NA, g.asr)) %>%
    #mutate(desc.rhodopsins = sapply(node, function(x) filter(phylophlan, node %in% getDescendants(phylophlan.phylo, x), !is.na(family)) %>% nrow)) %>%
    #mutate(g.asr = ifelse(desc.rhodopsins > 0, g.asr, NA)) %>%
    mutate(g = ifelse(!is.na(g.asr), g.asr, g)) %>%
    group_by(g, family) %>%
    arrange(-n()) %>%
    group_by(g) %>%
    mutate(dominant_family = first(na.omit(subfamily))) %>%
    mutate(g.asr = ifelse(is.na(dominant_family), NA, g.asr)) %>%
    `class<-`(c("tbl_tree", "data.frame")) %>%
    as.treedata

colors.habitats <- read.table(habitat_colors_file, sep = "\t", comment.char = "") %>%
    with(setNames(V2, V1))
colors.families <- read.table(family_colors_file, sep = "\t", comment.char = "") %>%
    with(setNames(V2, V1))

p <- ggtree(phylophlan.tree, aes(color = type), size = 0.2) +
    scale_color_manual(values = colors.habitats, na.value = "black") +
    new_scale_color() + new_scale_fill() +
    geom_treescale(width = 0.1) +
    geom_point2(aes(subset = ufboot >= 95 & SH_aLRT >= 80), shape = 15, color = "darkgray") +
    scale_size_continuous(range = c(1,2)) +
    # geom_text_repel(aes(label = cultured_organism), size = 3, hjust = 1, direction = "x", max.overlaps = Inf, min.segment.length = 0, seed = 42, box.padding = 0.5, segment.size = 0.2) +
    geom_tippoint(aes(subset = !is.na(subfamily), color = subfamily)) +
    geom_tiplab(aes(label = paste(fg, label), subset = isTip), align = T, size = 1, offset = 0.1, linesize = 0) +
    geom_tiplab(aes(label = cultured_organism, subset = !is.na(cultured_organism)), size = 3, align = T) +
    geom_text_repel(aes(label = f.asr), min.segment.length = 0, fontface = "bold") +
    geom_text_repel(aes(label = g.asr, color = dominant_family), min.segment.length = 0) +
    scale_colour_manual(values = colors.families)
    #geom_cladelab(mapping = aes(subset = rank == "f", node = node, label = label), offset = 0.05)
ggsave(output_file, p, width = 8, height = 12)
