#!/usr/bin/env Rscript
library(ggplot2)
library(readr)


args = commandArgs(trailingOnly=TRUE)


df <- read_tsv(args[1])
df$cluster <- as.factor(df$cluster)
size <- 0.2


p <- ggplot() + 
    geom_point(
        data=subset(df, cluster==-1),  # -1
        aes(x=c1, y=c2), size=size, color='grey90') +
    geom_point(
        data=subset(df, cluster!=-1), 
        aes(x=c1, y=c2, color=cluster), size=size) + 
    # geom_point(
    #     data=subset(anno, cluster==COI),  # -1
    #     aes(x=c1, y=c2), size=2*size, color='red') +
    guides(color=FALSE) +
    theme_classic() +
    # scale_x_continuous(limits=c(10, 25)) +
    # scale_y_continuous(limits=c(-5, 2)) +
    # scale_color_brewer(palette='Set2') +
    theme(axis.line=element_blank())


ggsave(args[2], p, width=20, height=20, units='cm')
