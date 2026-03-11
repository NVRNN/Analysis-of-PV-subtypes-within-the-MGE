library(tidyverse)
library(Seurat)

project_name <- "121023NathanielNapoli"


dat <- Read10X(data.dir = "/media/npnapoli/TOSHIBA EXT/Bioinformatics_Nestin/data/", unique.features = TRUE)
dat

grep(pattern = "PCDHA9", x = rownames(dat), value = TRUE)

dat %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "gene") %>%
  filter(grepl(pattern = "PCDHA9", x = gene)) %>%
  pivot_longer(-gene, names_to = "cell_barcode", values_to = "count") %>%
  group_by(gene) %>%
  summarise(total_counts = sum(count))
# Let's sum the counts for each gene



counts <- 
  dat %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "gene") %>%
  mutate(gene = str_match(string = gene, pattern = "^([^\\.]+)")[,2]) %>%
  pivot_longer(-gene, names_to = "cell_barcode", values_to = "count") %>%
  group_by(gene, cell_barcode) %>%
  summarise(count = sum(count)) %>%
  ungroup() %>%
  pivot_wider(names_from = cell_barcode, values_from = count)

# Save filtered data after regressing cell cycle
save(filename = "filtered.cell_cycle.h5Seurat",
             overwrite = TRUE)

