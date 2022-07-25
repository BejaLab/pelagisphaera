library(dplyr)
library(tidyr)
library(readxl)
library(tibble)
library(corrplot)
library(igraph)
library(rnaturalearth)
library(ggrepel)
library(ggnewscale)
library(tools)

with(snakemake@input, {
    salazar_file <<- salazar
    mgs_file     <<- mgs
    om_rgc_file  <<- om_rgc
    pelagi_file  <<- pelagi
    profile_file <<- profile
    usearch_file <<- usearch
    taxonomy_file <<- taxonomy
})
with(snakemake@output, {
    cogcor_file   <<- cogcor
    outliers_file <<- outliers
    pr_blh_file   <<- pr_blh
    biogeography_file <<- biogeography
    layers_file   <<- layers
})
with(snakemake@params, {
    Genus <<- genus
    gtdb_genus <<- gtdb_genus
})

taxonomy <- read.table(taxonomy_file, col.names = c("label", "taxon"), sep = "\t") %>%
    extract(taxon, into = "Order.gtdb", regex = "o__([^;]+)", remove = F) %>%
    extract(taxon, into = "Family.gtdb", regex = "f__([^;]+)", remove = F) %>%
    extract(taxon, into = "Genus.gtdb", regex = "g__([^;]+)") %>%
    extract(label, into = "Database", regex = "(GC.)_", remove = F)

samples <- read_excel(salazar_file, sheet = "Table_W1", .name_repair = "universal") %>%
    filter(MetaG.MetaT == "MetaG") %>%
    mutate(Layer = factor(Layer, levels = c("SRF", "DCM", "MES", "MIX")))
stations <- distinct(samples, Station, Longitude, Latitude)
stations_upper <- filter(samples, Layer %in% c("SRF", "DCM")) %>%
    distinct(Station, Longitude, Latitude)
coords <- select(samples, ENA_ID, ENA_Run_ID, Longitude, Latitude) %>%
    separate_rows(ENA_Run_ID, sep = "\\|") %>%
    {bind_rows(select(., sra = ENA_ID, Longitude, Latitude), select(., sra = ENA_Run_ID, Longitude, Latitude))}

Pelagi <- read.table(pelagi_file, header = T, fill = T, sep = "\t", na.string = "") %>%
    separate_rows(sra, sep = ", ") %>%
    left_join(coords, by = "sra") %>%
    mutate(Longitude = ifelse(is.na(Longitude), lon, Longitude), Latitude = ifelse(is.na(Latitude), lat, Latitude)) %>%
    # filter(!is.na(Latitude)) %>%
    mutate(Lat_rnd = round(Latitude), Lon_rnd = round(Longitude)) %>%
    distinct(label, organism, Lat_rnd, Lon_rnd, .keep_all = T) %>%
    group_by(organism)

target_labels <- filter(Pelagi, genus == Genus) %>%
    pull(label)

blast.names <- c("label", "query", "target", "identity", "length", "mismatches", "gaps", "qstart", "qend", "tstart", "tend", "evalue", "bitscore")
target_ids <- read.table(usearch_file, col.names = blast.names) %>%
    left_join(taxonomy, by = "label") %>%
    filter(!is.na(Order.gtdb)) %>%
    arrange(query, -identity) %>%
    distinct(query, .keep_all = T) %>%
    filter(Genus.gtdb %in% gtdb_genus, identity >= 95) %>%
    pull(query)

mgs <- read.table(mgs_file, col.names = c("COG", "description"), sep = "\t")

PR.COG <- filter(mgs, grepl("^PR", description)) %>%
    pull(COG)
blh.COG <- filter(mgs, grepl("^blh", description)) %>%
    pull(COG)

reads   <- read_excel(salazar_file, sheet = "Table_W2", .name_repair = "universal")

atlas <- read.csv(om_rgc_file) %>%
    rename(query = OM.RGC_ID)
COGs <- select(atlas, OMRGC_ID = query, OG)

Pelagi.poly <- filter(Pelagi, !is.na(Longitude)) %>%
    slice(chull(Longitude, Latitude))
Pelagi.labels <- filter(Pelagi, !is.na(Longitude)) %>%
    summarize(Longitude = mean(Longitude), Latitude = mean(Latitude), PR = first(PR))

profiles.raw <- read.csv(profile_file)
profiles <- gather(profiles.raw, PANGAEA.sample.id, profile, -OMRGC_ID) %>%
    left_join(COGs, by = "OMRGC_ID") %>%
    left_join(select(samples, PANGAEA.sample.id, Station, Layer), by = "PANGAEA.sample.id")

profiles_pelagi <- filter(profiles, OMRGC_ID %in% target_ids, Layer %in% c("DCM", "SRF")) %>%
    group_by(Station, OG) %>%
    summarize(profile = sum(profile), .groups = "drop") %>%
    spread(OG, profile) %>%
    column_to_rownames("Station")
profiles_layers <- filter(profiles, OMRGC_ID %in% target_ids) %>%
    group_by(PANGAEA.sample.id, Layer, Station, OG) %>%
    summarize(profile = sum(profile), .groups = "drop")

profiles.cor <- cor(profiles_pelagi)
pdf(cogcor_file, width = 10, height = 10)
corrplot(profiles.cor, order = "hclust")
dev.off()

graph <- profiles.cor %>% `[`(which(! rownames(.) %in% c(PR.COG, blh.COG)), which(! colnames(.) %in% c(PR.COG, blh.COG))) %>%
    graph.adjacency(weighted = T, mode = "lower")
comps <- delete.edges(graph, E(graph)[ weight < 0.9 ]) %>% decompose
biggest.group <- lapply(comps, vcount) %>%
    unlist %>%
    {which(. == max(.))}
cogs <- V(comps[[biggest.group]]) %>% names

outliers <- select(profiles_pelagi, !!cogs) %>%
    t %>% as.data.frame %>%
    mutate_all(function(x) (x - median(x))/mad(x)) %>%
    t %>% as.data.frame %>%
    gather(OG, outl) %>%
    filter(!is.nan(outl), !is.infinite(outl)) %>%
    group_by(OG) %>%
    mutate(outl.med = max(abs(outl))) %>%
    ungroup %>%
    arrange(-outl.med) %>%
    mutate(OG = factor(OG, levels = unique(OG)))
p <- ggplot(outliers, aes(x = OG, y = outl)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(outliers_file, p, width = 16, height = 6)

cogs.best <- distinct(outliers, OG) %>% tail(n = 20) %>% pull %>% as.character

profiles.all <- left_join(profiles.raw, COGs, by = "OMRGC_ID") %>%
    filter(OG %in% cogs.best) %>%
    group_by(OG) %>%
    summarize_at(vars(starts_with("TARA")), sum) %>%
    gather(PANGAEA.sample.id, profile, -OG) %>%
    left_join(samples, by = "PANGAEA.sample.id")

profiles_upper_stations <- filter(profiles.all, Layer %in% c("DCM", "SRF")) %>%
    group_by(Station) %>%
    summarize(MG.all = mean(profile, na.rm = T), .groups = "drop")
profiles_all_samples <- group_by(profiles.all, PANGAEA.sample.id) %>%
    summarize(MG.all = mean(profile, na.rm = T), .groups = "drop")

distinct_stations <- distinct(stations, Station, .keep_all = T)

MG.threshold <- 0.9
data <- select(profiles_pelagi, !!cogs.best, PR = !!PR.COG, blh = !!blh.COG) %>%
    mutate(MG = rowMeans(.[!!cogs.best]), MG.0 = rowSums(.[!!cogs.best] > 0)) %>%
    select(MG, PR, blh, MG.0) %>%
    mutate(operon = (PR + blh) / 2, MG.enough = MG.0 >= MG.threshold * length(cogs.best)) %>%
    mutate(operon.ratio = operon / MG, operon.ratio = ifelse(MG.enough, operon.ratio, ifelse(!is.nan(operon) & operon > 0, 1, 0))) %>%
    mutate(PR.ratio     =     PR / MG,     PR.ratio = ifelse(MG.enough, PR.ratio,     ifelse(!is.nan(PR)     & PR > 0,     1, 0))) %>%
    mutate(blh.ratio    =    blh / MG,    blh.ratio = ifelse(MG.enough, blh.ratio,    ifelse(!is.nan(blh)    & blh > 0,    1, 0))) %>%
    rownames_to_column("Station") %>%
    mutate(Station = as.numeric(Station)) %>%
    right_join(distinct_stations, by = "Station") %>%
    left_join(profiles_upper_stations, by = "Station")

p <- filter(data, blh.ratio >= 0, PR.ratio >= 0, blh.ratio > 0 | PR.ratio > 0, MG.enough) %>%
    ggplot(aes(x = blh, y = PR, color = log10(MG))) +
    geom_point() +
    coord_fixed() +
    scale_x_log10() +
    scale_y_log10() +
    geom_smooth() +
    theme_bw()
ggsave(pr_blh_file, p, width = 5, height = 5)

PR.colors <- c(yes = "darkblue", blh_only = "blue", no = "black", pseudo = "darkgray")

PR.col  <- "blue"
blh.col <- "#ffb380ff"

world <- ne_countries(scale = "medium", returnclass = "sf")

empty_stations <- left_join(stations_upper, data, by = "Station", suffix = c("", "y")) %>%
    filter(MG.0 == 0)
p <- ggplot(world) + geom_sf(color = "#ccffaa", fill = "#ccffaa") + coord_sf(expand = F) +
    geom_point(empty_stations, mapping = aes(y = Latitude, x = Longitude), shape = 3, size = 1) +
    geom_point(filter(data, MG.enough), mapping = aes(y = Latitude, x = Longitude, color = operon.ratio, size = log10(MG / MG.all)), shape = 16) +
    scale_colour_gradient2(high = "darkblue", mid = "gray", low = "black") + new_scale_color() +
    geom_point(filter(data, !MG.enough), mapping = aes(y = Latitude, x = Longitude, color = operon.ratio, size = log10(MG / MG.all)), shape = 1) +
    scale_colour_gradient2(high = "darkblue", mid = "gray", low = "black") + new_scale_color() +
    scale_size(range = c(1, 8)) +
    geom_label_repel(Pelagi.labels, mapping = aes(y = Latitude, x = Longitude, label = organism, color = PR),
        arrow = arrow(length = unit(0.005, "npc")), box.padding = 1, min.segment.length = 0, max.overlaps = Inf, fill = NA) +
    scale_color_manual(values = PR.colors) +
    geom_polygon(Pelagi.poly, mapping = aes(y = Latitude, x = Longitude, group = organism), alpha = 0.1, color = "black") +
    theme_bw()
ggsave(biogeography_file, p, device = cairo_pdf, width = 12, height = 7)

dbl_samples <- group_by(samples, Station) %>%
    filter(Layer != "MIX") %>%
    filter(n_distinct(Layer) > 1) %>%
    pull(PANGAEA.sample.id)

# NB: a proper solution is to sum over Station+Layer

rm.sample <- "TARA_B100000959"
pelagi_data <- filter(profiles_layers) %>%
    mutate(group = case_when(OG %in% cogs.best ~ "MGs", OG %in% c(PR.COG, blh.COG) ~ "operon")) %>%
    # filter(PANGAEA.sample.id != rm.sample) %>%
    group_by(PANGAEA.sample.id) %>%
    mutate(MG.0 = sum(group == "MGs" & profile > 0, na.rm = T)) %>%
    group_by(group, PANGAEA.sample.id) %>%
    summarize(profile = mean(profile), MG.0 = first(MG.0), .groups = "drop") %>%
    spread(group, profile) %>%
    select(-`<NA>`) %>%
    left_join(profiles_all_samples, by = "PANGAEA.sample.id") %>%
    left_join(samples, by = "PANGAEA.sample.id") %>%
    filter(Layer != "MIX") %>%
    group_by(Station, Layer) %>%
    mutate(Layer = factor(Layer, levels = c("SRF", "DCM", "MES")), OS.region = factor(OS.region), Station = factor(Station)) %>%
    mutate(MG.enough = MG.0 >= MG.threshold * length(cogs.best)) %>%
    mutate(MG.ratio = MGs / MG.all, operon.ratio = operon / MGs, operon.ratio = ifelse(MG.enough, operon.ratio, ifelse(!is.nan(operon) & operon > 0, 1, 0)))

p <- ggplot(filter(pelagi_data, MGs > 0), aes(x = Layer, y = MG.ratio)) +
    geom_violin(trim = F, scale = "count") +
    geom_jitter(aes(color = operon.ratio, shape = MG.enough), width = 0.05, size = 5) +
    scale_colour_gradient2(high = "darkblue", low = "black", mid = "gray") +
    scale_y_log10() +
    scale_shape_manual(values = c(1, 16)) +
    ylab("Ca. Pelagisphaera fraction") + xlab("Layer") +
    theme_bw()
ggsave(layers_file, p, device = cairo_pdf, width = 6, height = 6)
