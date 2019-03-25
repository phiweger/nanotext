#!/usr/bin/env Rscript
library(ggplot2)
library(readr)


args = commandArgs(trailingOnly=TRUE)


# Clostridia
# Gammaproteobacteria
# Chlamydiia
# Oxyphotobacteria
anno <- read_tsv(args[1])
size <- 0.5
p <- ggplot() + 
    geom_point(
        data=subset(anno, is.na(anno$subclade)),  # -1
        aes(x=c1, y=c2), size=size, color='grey90') +
    geom_point(
        data=subset(anno, !is.na(anno$subclade)), 
        aes(x=c1, y=c2, color=subclade), size=size) + 
    # geom_point(
    #     data=subset(anno, cluster==COI),  # -1
    #     aes(x=c1, y=c2), size=2*size, color='red') +
    # guides(color=FALSE) +
    theme_classic() +
    # scale_x_continuous(limits=c(10, 25)) +
    # scale_y_continuous(limits=c(-5, 2)) +
    scale_color_brewer(palette='Set2') +
    theme(axis.line=element_blank())
ggsave(paste0(args[2], 'ecotypes.subclade.pdf'), p, 
    height=8, width=15, units='cm')


anno <- read_tsv(args[1])
size <- 0.5
r <- ggplot() + 
    geom_point(
        data=subset(anno, is.na(anno$clade)),  # -1
        aes(x=c1, y=c2), size=size, color='grey90') +
    geom_point(
        data=subset(anno, !is.na(anno$clade)), 
        aes(x=c1, y=c2, color=clade), size=size) + 
    # geom_point(
    #     data=subset(anno, cluster==COI),  # -1
    #     aes(x=c1, y=c2), size=2*size, color='red') +
    # guides(color=FALSE) +
    theme_classic() +
    # scale_x_continuous(limits=c(10, 25)) +
    # scale_y_continuous(limits=c(-5, 2)) +
    scale_color_brewer(palette='Set2') +
    theme(axis.line=element_blank())
ggsave(paste0(args[2], 'ecotypes.clade.pdf'), r, 
    height=8, width=15, units='cm')


anno <- read_tsv(args[1])
anno$cluster <- as.factor(anno$cluster)
size <- 1
s <- ggplot() + 
    geom_point(
        data=subset(anno, cluster==-1),  # -1
        aes(x=c1, y=c2), size=size, color='grey90') +
    geom_point(
        data=subset(anno, cluster!=-1), 
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
ggsave(paste0(args[2], 'ecotypes.clusters.pdf'), s, 
    height=8, width=15, units='cm')


# size <- 0.5
# q <- ggplot() + 
#     geom_point(
#         data=subset(anno, is.na(anno$clade)),  # -1
#         aes(x=c1, y=c2), size=size, color='grey90') +
#     geom_point(
#         data=subset(anno, !is.na(anno$clade)), 
#         aes(x=c1, y=c2, color=clade), size=size) + 
#     # geom_point(
#     #     data=subset(anno, cluster==COI),  # -1
#     #     aes(x=c1, y=c2), size=2*size, color='red') +
#     # guides(color=FALSE) +
#     facet_wrap(~collection) +
#     theme_classic() +
#     theme(strip.background=element_blank()) +
#     scale_x_continuous(limits=c(10, 25)) +
#     scale_y_continuous(limits=c(-5, 2)) +
#     scale_color_brewer(palette='Set2') +
#     theme(axis.line=element_blank())
# ggsave('ecotypes_model_i_clades_facet.pdf', q, height=6, width=22, units='cm')