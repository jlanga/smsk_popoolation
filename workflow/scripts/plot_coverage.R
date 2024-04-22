#!/usr/bin/env Rscript
library(tidyverse)
library(tidyquant)


frequencies <-
  list.files(
    path = "results/preprocess/coverage",
    pattern = ".hist$",
    full.names = TRUE
  ) %>%
  map(
    function(x) read_tsv(
      file = x,
      col_names = c("population", "coverage", "frequency"),
      col_types = "cii"
    )
  ) %>%
  bind_rows()


coverages <- frequencies %>%
  filter(coverage > 0) %>%
  ggplot(aes(x = coverage, y = frequency)) +
  geom_ma(ma_fun = SMA, linetype="solid", color = "black") +
  scale_y_log10() +
  facet_wrap(~population)
ggsave("results/preprocess/coverage/coverage.pdf", width = 297, height = 210, units = "mm")



frequencies %>%
  filter(coverage > 10) %>%
  group_by(population) %>%
  filter(frequency == max(frequency)) %>%
  mutate(
    max_coverage = 1.5 * coverage,
    min_coverage = 0.5 * coverage
  ) %>%
  write_tsv("results/preprocess/coverage/coverage.tsv")
