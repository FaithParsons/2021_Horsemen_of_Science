library(tidyverse)
library(ggplot2)
library(janitor)
library(ggthemes)


horses <- readxl::read_xlsx("../data/HR BEHAV database for Faith.xlsx", sheet="data") %>%
  filter(timepoint == 2 | timepoint==3 | timepoint == 4) %>%
  mutate(tp = recode(timepoint, "tp2" = 2, "tp3" = 3, "tp4" = 4)) %>%
  as.data.frame()

ggplot (horses, aes(x=factor(timepoint), y=SniffPC,fill=factor(stimulus))) +
  geom_boxplot() +
  geom_jitter() +
  theme_classic()

ggplot (horses, aes(x=factor(timepoint), y=ChewPC,fill=factor(stimulus))) +
  geom_boxplot() +
  geom_jitter() +
  theme_classic()

ggplot (horses, aes(x=factor(timepoint), y=EarsPC,fill=factor(stimulus))) +
  geom_boxplot() +
  geom_jitter() +
  theme_classic()



ggplot (horses, aes(x=factor(timepoint), y=WithdrPC,fill=factor(stimulus))) +
  geom_boxplot() +
  geom_jitter() +
  theme_classic()


ggplot (horses, aes(x=factor(timepoint), y=InteractPC,fill=factor(stimulus))) +
  geom_boxplot() +
  geom_jitter() +
  theme_classic()



