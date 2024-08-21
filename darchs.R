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

good <- good |>
  mutate(slength2 = standarize(length),
         deam = str_detect(arch_code, "å"))

good <- good |>
  group_by(deam) |>
  mutate(slength_deam = standarize(length))


good <- good |>
  select(-slength, -slength2) |>
  rename( slength = slength_deam)

view(good)

gen_type()

x <- archs_pid$arch
Bdeam <- str_detect(x, deam)
Blxg <- str_detect(x, lxg)
Bendo <- str_detect(x, endo)
Bdeam_lxg <- Bdeam & Blxg

type_pid <- case_when(
  Bdeam_lxg ~ "Deam_LXG",
  Bdeam ~ "deam",
  Bendo ~ "EndoV",
  .default = "."
)

archs_pid <- archs_pid |>
  mutate(type = type_pid)

archs_pid |>
  left_join(a)

blasts <- read_tsv("blasts.tsv")
names(blasts)
archs_pid <- blasts |>
  distinct(pid, bitscore) |>
  right_join(archs_pid, join_by(pid))


archs_pid <- archs_pid |>
  mutate(sbitscore = standarize(bitscore)) |>
  relocate(pid, type, arch_code, arch, length, slength, sbitscore, bitscore)

d_lv <- stringdistmatrix(archs_pid$arch_code, archs_pid$arch_code, method = "lv")
d_cos <- stringdistmatrix(archs_pid$arch_code, archs_pid$arch_code, method = "cosine")
d_jac <- stringdistmatrix(archs_pid$arch_code, archs_pid$arch_code, method = "jaccard")

hm <- function(m) {
  m[1:10, 1:16]
}

lev <- as_tibble(d_lv) |>
  mutate(across(everything(), as.integer))
names(lev) <- archs_pid$pid

lev <- lev |>
  mutate(
    pid = archs_pid$pid) |>
  relocate(pid)

slev <- lev |>
  mutate(across(!pid, standarize))

slev <- slev |>
  mutate(
    slegth = archs_pid$slength,
    sbitscore = archs_pid$sbitscore,
    type = archs_pid$type) |>
  relocate(pid, type)

# https://www.r-bloggers.com/2019/05/quick-and-easy-t-sne-analysis-in-r/
library(M3C)

sub_num <- slev |>
  map_lgl(is.numeric)
slev_numeric <- slev[sub_num]

tsne(slev_numeric, labels=as.factor(slev$type))

# https://datavizpyr.com/how-to-make-tsne-plot-in-r/
library(Rtsne)
theme_set(theme_bw(18))
set.seed(142)

tSNE_fit <- Rtsne(slev_numeric)

tSNE_df <- tSNE_fit$Y %>%
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(ID=row_number())

p <- tSNE_df %>%
  ggplot(aes(x = tSNE1,
             y = tSNE2,
             color = type
  geom_point(alpha = 1/3)+
  theme(legend.position="bottom")+
  title("tSNE Bacillota Blast YwqJ, YwqL")
p
ggsave(p, "tSNE_proteins_bacillota.svg")
