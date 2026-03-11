library(tidyverse)

project_name <- "230622NathanielNapoli"
data_dir <- paste("~/localwork/facility/data", project_name, sep = "/")
out_dir <- paste("~/localwork/facility/analysis_output", project_name, sep = "/")
dir.create(out_dir)

meta <- read_tsv(file = paste(data_dir, "GSE109796_series_matrix.txt", sep = "/"), skip = 27)


# Fix column names
meta %>% colnames()

meta <-
  meta %>%
  rename_with(~ gsub(pattern = "-", replacement = "_", x = .)) %>%
  rename(data_type = 1)

meta


meta <-
  meta %>%
  pivot_longer(cols = -data_type, names_to = "sample") %>%
  filter(grepl(pattern = "age:", x = value) | grepl(pattern = "tissue:", x = value)) %>%
  filter(!grepl(pattern = "Brain", x = value)) %>%
  select(-data_type) %>%
  mutate(info_type = str_split_fixed(string = value, pattern = ":", n = 2)[,1],
         value = str_split_fixed(string = value, pattern = " ", n = 2)[,2]) %>%
  select(sample, info_type, value) %>%
  pivot_wider(names_from = info_type, values_from = value)

meta

meta %>% 
  select(tissue, age) %>%
  table()
# Do these values match up to the paper?


# Try another approach - is this the "correct" way?
# https://bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html
library(GEOquery)

test <- getGEO("GSE109796") # Try GSEMatrix = FALSE
class(test)
length(test)

class(test$GSE109796_series_matrix.txt.gz)
test$GSE109796_series_matrix.txt.gz$title

test2 <- tibble(sample = test$GSE109796_series_matrix.txt.gz$title,
                tissue = test$GSE109796_series_matrix.txt.gz$characteristics_ch1.1,
                age = test$GSE109796_series_matrix.txt.gz$characteristics_ch1.2)
test2

test2 <-
  test2 %>%
  mutate(tissue = str_split_fixed(pattern = " ", string = tissue, n = 2)[,2],
         age = str_split_fixed(pattern = " ", string = age, n = 2)[,2])
test2
meta

all_equal(meta, test2)
class(test2)
class(meta)

head(test2)
head(meta)


test2 %>%
  write_tsv(file = paste(out_dir, "metadata.tsv", sep = "/"))
