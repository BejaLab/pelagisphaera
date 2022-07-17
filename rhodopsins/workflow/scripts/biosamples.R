library(dplyr)
library(tidyr)
library(XML)

xml_files   <- unlist(snakemake@input)
ext         <- unlist(snakemake@params["ext"])
output_file <- unlist(snakemake@output)
lapply(xml_files, xmlParse) %>%
    lapply(xmlToList) %>%
    lapply(unlist) %>%
    setNames(gsub(ext, "", basename(xml_files))) %>%
    bind_rows(.id = "label") %>%
    unite(key,   starts_with("DocumentSummary.SampleData.BioSample.Attributes.Attribute..attrs.attribute_name"), sep = "@@@") %>%
    unite(value, starts_with("DocumentSummary.SampleData.BioSample.Attributes.Attribute.text"),          sep = "@@@") %>%
    rename(OrganismName = DocumentSummary.SampleData.BioSample.Description.Organism.OrganismName) %>%
    rename(taxonomy_id  = DocumentSummary.SampleData.BioSample.Description.Organism..attrs.taxonomy_id) %>%
    rename(taxonomy_name = DocumentSummary.SampleData.BioSample.Description.Organism..attrs.taxonomy_name) %>%
    select(label, Organism = DocumentSummary.Organism, taxonomy_id, taxonomy_name, key, value) %>%
    separate_rows(key, value, sep = "@@@") %>%
    mutate(key = gsub("[ -]", "_", tolower(key)), value = tolower(value)) %>%
    filter(! value %in% c("na", "not applicable", "missing", "false", "not provided", "not collected")) %>%
    mutate(is.metagenome  = grepl("metagenom", value) | grepl("metagenom", key) | grepl("uncultured", value) | grepl("\\b(SCGC)\\b", Organism)) %>%
    mutate(is.marine      = grepl("(marine|seawater|ocean|sea water|sea beach|arenae)", value)) %>%
    mutate(is.soil    = grepl("\\b(soil|peat|mud|steppe|terrae|forest)\\b", value)) %>%
    mutate(is.salt_lake   = grepl("\\b(soda lake|hypersaline|salt lake)\\b", value)) %>%
    mutate(is.freshwater  = grepl("\\b(freshwater|lentic|river|pond|groundwater|springs|spring|lake)\\b", value)) %>%
    mutate(is.from_animal = grepl("\\b(shelfordella|xiphinema|mouse|python|reticulitermes|cephalotes|gallus|pig|termite|bos|sus|goat|sheep|deer|gut|bodily|fecal|cattle|feces|faeces|faecal|stool|excreta|yak|homo|gastrointestinal|ruminant|buffalo|human|mus musculus|caecum)(-associated)?\\b", value)) %>%
    mutate(is.from_marine_inv = grepl("\\b(lissoclinum|polychaeta|petrosia|vazella)\\b", value)) %>%
    mutate(is.technogen   = grepl("\\b(bioreactor|sludge|wastewater|digester)\\b", value)) %>%
    mutate(is.symbiont    = grepl("^host", key) | grepl("\\b(symbiont|symbiotic)\\b", value)) %>%
    mutate(is.water       = grepl("\\b(aquatic|water|aquifer)\\b", value)) %>%
    mutate(is.cultured    = grepl("culture_collection", key) | grepl("\\b(culture|type strain)\\b", value) | !grepl("\\bbacterium\\b", Organism) & key == "strain") %>%
    group_by(label, Organism) %>%
    summarize(type = case_when(
            any(is.from_animal,     na.rm = T) ~ "terrestrial animals",
            any(is.from_marine_inv, na.rm = T) ~ "marine invertebrates",
            any(is.soil,            na.rm = T) ~ "soil",
            any(is.salt_lake,       na.rm = T) ~ "salt lake",
            any(is.freshwater,      na.rm = T) ~ "freshwater",
            any(is.marine,          na.rm = T) ~ "marine"
        ),
        is.metagenome = any(is.metagenome, na.rm = T),
        is.cultured = !is.metagenome & any(is.cultured, na.rm = T)
    ) %>%
    write.table(file = output_file, col.names = T, row.names = F)
