#!/usr/bin/env Rscript

library(tidyverse)
library(stringr)
library(stringdist)
# library(GenomicRanges)

l <- length

IN <- "iscan.tsv"

iscan <- read_tsv(IN)

archs <- iscan |>
  filter(analysis == "Pfam", recommended) |>
  group_by(pid) |>
  reframe(
    domain = memberDB, start = start, end = end,
    length = length, domain_txt = memberDB_txt
  ) |>
  arrange(pid, start, end) |>
  mutate(
    start = as.integer(start),
    end = as.integer(end),
    length = as.integer(length)
  )

# Short Archs
Sarchs <- archs |>
  select(-length, -domain_txt)


# List Archs
Larchs <- Sarchs |>
  group_by(pid) |>
  group_split()

names(Larchs) <- Larchs |>
  map_chr(\(x) unique(x$pid))

# Asserts
assert_increasing_start <- Larchs |>
  map_lgl(\(x) reduce(c(0, x$start), `<`)) |>
  all()

assert_increasing_end <- Larchs |>
  map_lgl(\(x) reduce(c(0, x$end), `<`)) |>
  all()

stopifnot(
  "non-strictly increasing starts" =
    assert_increasing_start,
  "non-strictly increasing ends" =
    assert_increasing_end
)

archs_pid <- Larchs |>
  map(\(x) tibble(
    pid = unique(x$pid),
    arch = str_flatten(x$domain, collapse = ",")
  )) |>
  do.call(bind_rows, args = _)

rnorm(n = 10)

standarize <- function(x) {
  mean_x <- mean(x)
  sd_x <- sd(x)
  (x - mean_x) / sd_x
}





one_lettercode <- function(doms) {
  doms <- unique(doms)

  OG <- c(46, 60:70, 97:122)
  START <- 192
  OFF <- 64 # para hacerlos todavia mas distintos

  if (length(OG) >= length(doms)) {
    SAMPLE <- sample(OG, length(doms), replace = F)
  } else {
    extra <- START:(START + (length(doms) - length(OG)) + OFF)
    SAMPLE <- sample(c(OG, extra), length(doms), replace = F)
  }

  START_U <- 192 + length(doms)

  STEP <- #
    OUT <- vector(mode = "list", length = length(doms))
  names(OUT) <- doms

  i <- 1
  for (dom in doms) {
    Ucode <- intToUtf8(SAMPLE[[i]])
    OUT[[dom]] <- Ucode
    i <- i + 1
  }
  OUT
}

ONE_LETTER <- one_lettercode(archs$domain)

replace_to_oneletter <- function(archs, code) {
  keys <- names(code)
  for (key in keys) {
    archs <- str_replace_all(archs, key, code[[key]])
  }
  str_replace_all(archs, ",", "")
}


archs_pid <- archs |>
  distinct(pid, length) |>
  right_join(archs_pid) |>
  mutate(slength = standarize(length))

archs_pid <- archs_pid |>
  mutate(arch_code = replace_to_oneletter(arch, ONE_LETTER))

lxg <- "PF04740"
deam <- "PF14431"
endo <- "PF04493"


bad <- archs_pid |>
  filter(!str_detect(arch, endo) & !str_detect(arch, deam))

good <- archs_pid |>
  filter(str_detect(arch, endo) |
    (str_detect(arch, lxg) & str_detect(arch, deam)))
