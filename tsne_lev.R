library(tidyverse)
library(ggthemes)
set.seed(142)

# https://www.r-bloggers.com/2019/05/quick-and-easy-t-sne-analysis-in-r/
# library(M3C)

# https://datavizpyr.com/how-to-make-tsne-plot-in-r/
library(Rtsne)
theme_set(theme_bw(18))
l <- length

# Already normalized
IN <- "slev.tsv"
slev <- read_tsv(IN)

# sub_num <- slev |>
#   map_lgl(is.numeric)
# slev_numeric <- slev[sub_num]

r_slev <- seq_len(nrow(slev))
r_sampled <- sample(r_slev, 2048, replace = F)
perplexities <- c(5, 30, 50)


INPUT <- slev[r_sampled, 3:l(slev)]
#
# tSNE_fit <- Rtsne(slev[3:l(slev)],
#   check_duplicates = FALSE,
#   perplexity = 156,
#   num_threads = 0,
#   normalize = F
# )

tSNE_df <- tSNE_fit$Y %>%
  as.data.frame() %>%
  rename(
    tSNE1 = "V1",
    tSNE2 = "V2"
  ) %>%
  mutate(
    pid = slev$pid,
    type = slev$type
  )

# write_tsv(tSNE_df, "156p.tsne.tsv")



p <- tSNE_df %>%
  ggplot(aes(
    x = tSNE1,
    y = tSNE2,
    color = type
  )) +
  geom_point(alpha = 1 / 4, size = 4, shape = 1) +
  labs(
    title = "YwqJ & YwqL Domain Architecture Groups",
    subtitle = "Blast against all Bacillota phylum\nt-SNE using Levenshtein distance",
    caption = "author: Becerra-Soto E."
  ) +
  xlab("tSNEx") +
  ylab("tSNEy") +
  scale_color_hue(labels = c("Other", "Deaminase", "Deaminase & LXG", "Endonuclease V")) +
  theme_fivethirtyeight() +
  coord_equal(ratio = 1)

p <- p + guides(color = guide_legend(
  title = "Group",
  legend.title = element_text(
    size = 16,
    face = "italic"
  ),
  override.aes = list(alpha = 1)
))

#
# ggsave("darchs_tSNE_bacillota.p156.svg", p, width = 12, height = 10)
