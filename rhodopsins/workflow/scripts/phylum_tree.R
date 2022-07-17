
library(dplyr)
library(tidyr)
library(ape)
library(ggtree)
library(igraph)
library(treeio)
library(ggnewscale)
library(tibble)
library(ggrepel)

with(snakemake@input, {
    biosamples_file <<- biosamples
    gtdbtk_tree     <<- gtdbtk_tree
    gtdbtk_taxonomy <<- gtdbtk_taxonomy
    rhodopsins_file <<- rhodopsins
    family_colors_file <<- family_colors
    habitat_colors_file <<- habitat_colors
})
with(snakemake@params, {
    shallow.threshold <<- shallow_threshold
    outlier.threshold <<- outlier_threshold
    phylum <<- phylum
})
output_file <- unlist(snakemake@output)

metadata <- read.table(biosamples_file, header = T)

rhodopsins <- read.table(rhodopsins_file, header = T) %>%
    group_by(family) %>%
    mutate(n = n()) %>%
    arrange(-n) %>%
    ungroup %>%
    distinct(label, .keep_all = T)

gtdbtk <- read.tree(gtdbtk_tree)
node.id <- with(gtdbtk, which(grepl(phylum, node.label)))
verruco.full <- extract.clade(gtdbtk, length(gtdbtk$tip.label) + node.id)

edges <- with(verruco.full, data.frame(edge, edge.length))
outlier.tips <- nodepath(verruco.full) %>%
    lapply(function(x) filter(edges, X1 %in% x, X2 %in% x)) %>%
    lapply(pull, "edge.length") %>%
    lapply(sum) %>%
    {data.frame(dist = unlist(.), label = verruco.full$tip.label)} %>%
    filter(dist > outlier.threshold) %>%
    pull(label)

tip.clust <- cophenetic(verruco.full) %>%
    `<`(shallow.threshold) %>%
    graph_from_adjacency_matrix %>%
    clusters %>%
    with(data.frame(label.orig = names(membership), cluster = membership)) %>%
    mutate(is.ref = grepl("^(GB|RS)_", label.orig)) %>%
    mutate(label = sub("^(GB|RS)_", "", label.orig)) %>%
    mutate(is.refseq = grepl("^GCF_", label)) %>%
    left_join(rhodopsins, by = "label") %>%
    left_join(metadata, by = "label") %>%
    filter(!is.na(Organism)) %>%
    arrange(!is.ref, !is.refseq, !is.cultured, cluster) %>%
    distinct(cluster, .keep_all = T)

taxonomy <- read.table(gtdbtk_taxonomy, col.names = c("label1", "label2", "classification"), header = T, sep = "\t")
taxonomy <- bind_rows(
         select(taxonomy, label = label1, classification),
         select(taxonomy, label = label2, classification)
    ) %>%
    distinct(label, classification) %>%
    separate_rows(classification, sep = ";") %>%
    separate(classification, into = c("rank","taxon"), sep = "__") %>%
    spread(rank, taxon)

verruco.tree <- with(verruco.full, tip.label[! tip.label %in% tip.clust$label.orig | tip.label %in% outlier.tips]) %>%
    {ape::drop.tip(phy = verruco.full, tip = .)}
verruco <- as_tibble(verruco.tree) %>%
    mutate(is.tip = label %in% verruco.tree$tip.label) %>%
    mutate(label = sub("^(GB|RS)_", "", label)) %>%
    left_join(tip.clust, by = "label") %>%
    left_join(taxonomy, by = "label") %>%
    separate(label, into = c("support", "label"), sep = ":", convert = T, fill = "left") %>%
    separate_rows(label, sep = ";") %>%
    distinct(node, .keep_all = T) %>%
    separate(label, into = c("rank", "label"), sep = "__", extra = "merge", fill = "left") %>%
    replace_na(list(rank = "")) %>%
    mutate(g = ifelse(rank == "g", label, g)) %>%
    group_by(o) %>%
    mutate(mrca = getMRCA(verruco.tree, node[is.tip]) %>% replace(is.null(.), NA), mrca = ifelse(is.na(mrca), node, mrca)) %>%
    group_by(mrca) %>%
    mutate(is.mrca = node == mrca & sum(is.tip) > 0, rank = ifelse(is.mrca, "o", rank), label = ifelse(is.mrca, first(na.omit(o)), label), o = ifelse(is.mrca, label, o)) %>%
    group_by(o) %>%
    mutate(c = first(na.omit(c))) %>%
    group_by(f) %>%
    mutate(mrca = getMRCA(verruco.tree, node[is.tip]) %>% replace(is.null(.), NA), mrca = ifelse(is.na(mrca), node, mrca)) %>%
    group_by(mrca) %>%
    mutate(is.mrca = node == mrca & sum(is.tip) > 0, rank = ifelse(is.mrca, "f", rank), label = ifelse(is.mrca, first(na.omit(f)), label), f = ifelse(is.mrca, label, f)) %>%
    group_by(f) %>%
    mutate(dominant_subfamily_f = table(subfamily) %>% sort(d = T) %>% names %>% first(d = NA), dominant_subfamily_f_label = ifelse(is.na(dominant_subfamily_f) | rank != "f", "", f)) %>%
    group_by(g) %>%
    mutate(dominant_subfamily_g = table(subfamily) %>% sort(d = T) %>% names %>% first(d = NA), dominant_subfamily_g_label = ifelse(is.na(dominant_subfamily_g) | rank != "g", "", g)) %>%
    ungroup %>%
    mutate(dominant_subfamily = case_when(rank == "f" ~ dominant_subfamily_f, rank == "g" ~ dominant_subfamily_g, T ~ NA_character_),
           dominant_subfamily_label = ifelse(is.na(dominant_subfamily) | ! rank %in% c("g", "f"), "", label)
    ) %>%
    `class<-`(c("tbl_tree", "data.frame"))

#verruco.heatmap <- filter(verruco, !is.na(Organism)) %>%
#    select(label, type) %>%
#    column_to_rownames("label")

colors.habitats <- read.table(habitat_colors_file, sep = "\t", comment.char = "") %>%
    with(setNames(V2, V1))
colors.families <- read.table(family_colors_file, sep = "\t", comment.char = "") %>%
    with(setNames(V2, V1))

p <- ggtree(as.treedata(verruco), aes(color = type), size = 0.1, layout = "circular") +
    scale_color_manual(values = colors.habitats, na.value = "black") +
    new_scale_color() + new_scale_fill() +
    geom_tiplab(mapping = aes(subset = is.cultured), color = "red", label = "✱", size = 2) +
    geom_tippoint (mapping = aes(subset = !is.na(subfamily), color = subfamily)) +
    geom_cladelab (mapping = aes(subset = rank == "c", node = node, label = label, angle = angle), offset = 0.05) +
    scale_color_manual(values = colors.families, na.value = "gray") +
    geom_cladelab (mapping = aes(subset = rank == "o" & c == "Verrucomicrobiae", node = node, label = label)) +
    geom_highlight(mapping = aes(subset = rank == "o", fill  = label), alpha = 0.1) +
    geom_text_repel(aes(label = dominant_subfamily_label, color = dominant_subfamily, fontface = recode(rank, g = "plain", f = "bold")), max.overlaps = 100, segment.linetype = 2, segment.size = 0.2, min.segment.length = 0, point.padding = 0, force = 10) +
    geom_treescale(width = 0.2)
ggsave(output_file, p, width = 12, height = 12)
