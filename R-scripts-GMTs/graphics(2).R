pPpi <- ggplot(ppiV, aes(x = PPI, y = AUC, 
                         fill = PPI)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  labs(title = "Protein-Protein Interactions: AUC with ARCHS4 Human", y = "AUC") +
  theme_classic() + 
  # theme(plot.title = element_text(size = 20, hjust = 0.5), 
  #       axis.title = element_text(size = 16),
  #       axis.text = element_text(size = 14)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 11),
        plot.title = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "none",
        plot.margin = margin(2, 2, 2, 2, "cm"),
        plot.background = element_rect(
          fill = "white",
          colour = "black",
          size = 1)
        ) +
  scale_fill_brewer(palette = "OrRd")


# function to get the legend - got from www.sthda.com 
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# violin plot for ChEA

pCh <- ggplot(cheaV, aes(x = Expression, y = AUC, 
                                fill = Expression)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  labs(title = "ChEA 2016", y = "AUC") +
  theme_classic() + 
  # theme(plot.title = element_text(size = 20, hjust = 0.5), 
  #       axis.title = element_text(size = 16),
  #       axis.text = element_text(size = 14)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "top",
        legend.title = element_blank()) +
  # scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_fill_brewer(palette = "Paired", labels = c("CCLE", "GTEX", "ARCHS4 Human", "ARCHS4 Mouse"))

# violin plot for ENCODE

pEnc <- ggplot(encV, aes(x = Expression, y = AUC, 
                                   fill = Expression)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  labs(title = "ENCODE TF ChIP-seq 2015", y = "AUC") +
  theme_classic() + 
  # theme(plot.title = element_text(size = 20, hjust = 0.5), 
  #       axis.title = element_text(size = 16),
  #       axis.text = element_text(size = 14)) + 
  # scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position="none", 
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = 10))

pGBio <- ggplot(goBioV, aes(x = Expression, y = AUC, 
                                   fill = Expression)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  labs(title = "GO Biological Process 2017", y = "AUC") +
  theme_classic() + 
  # theme(plot.title = element_text(size = 20, hjust = 0.5), 
  #       axis.title = element_text(size = 16),
  #       axis.text = element_text(size = 14)) + 
  # scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position="none", 
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = 10))

pGCell <- ggplot(goCellV, aes(x = Expression, y = AUC, 
                                     fill = Expression)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  labs(title = "GO Cellular Component 2017", y = "AUC") +
  theme_classic() + 
  # theme(plot.title = element_text(size = 20, hjust = 0.5), 
  #       axis.title = element_text(size = 16),
  #       axis.text = element_text(size = 14)) + 
  # scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position="none", 
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = 10))

pGMol <- ggplot(goMolV, aes(x = Expression, y = AUC, 
                                   fill = Expression)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  labs(title = "GO Molecular Function 2017", y = "AUC") +
  theme_classic() + 
  # theme(plot.title = element_text(size = 20, hjust = 0.5), 
  #       axis.title = element_text(size = 16),
  #       axis.text = element_text(size = 14)) + 
  # scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position="none", 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = 10))

pKea <- ggplot(keaV, aes(x = Expression, y = AUC, 
                                fill = Expression)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  labs(title = "KEA 2015", y = "AUC") +
  theme_classic() + 
  # theme(plot.title = element_text(size = 20, hjust = 0.5), 
  #       axis.title = element_text(size = 16),
  #       axis.text = element_text(size = 14)) + 
  # scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position="none", 
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = 10))

pKgg <- ggplot(keggV, aes(x = Expression, y = AUC, 
                                 fill = Expression)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  labs(title = "KEGG 2016", y = "AUC") +
  theme_classic() + 
  # theme(plot.title = element_text(size = 20, hjust = 0.5), 
  #       axis.title = element_text(size = 16),
  #       axis.text = element_text(size = 14)) + 
  # scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position="none", 
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = 10))

pMgi <- ggplot(mgiV, aes(x = Expression, y = AUC, 
                                fill = Expression)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  labs(title = "MGI Mammalian Phenotype Lvl 4", y = "AUC") +
  theme_classic() + 
  # theme(plot.title = element_text(size = 20, hjust = 0.5), 
  #       axis.title = element_text(size = 16),
  #       axis.text = element_text(size = 14)) + 
  # scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position="none", 
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = 10))

# Use the code below to plot the violin plots
# Before doing so, get the legend from the pCh plot using the get_legend function at the top of the script
# Remove the legend from the pCh plot using pCh = pCh + theme(legend.position = "none")
# Then use the code below to plot with everything nicely arranged
# grid.arrange(pCh, pEnc, pGBio, pGCell, pGMol, pKea, pKgg, pMgi,legend,
#              ncol = 4, nrow = 3, layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 9, 9, 9)),
#              heights = c(2.5, 2.5, .5))

# This made a violin plot for the protein-protein interactions - run entirely separately.
# pProt = ggplot(prot_violinData, aes(x = PPI, y  = AUC, fill = PPI)) + 
#   geom_violin(trim = T, alpha = 0.5) + 
#   ggtitle("Protein-Protein Interactions") + 
#   geom_boxplot(outlier.shape = NA,width=0.1) + 
#   scale_y_continuous(limits = c(0.2,1.0)) +
#   theme(plot.title = element_text(size = 16),
#         axis.title = element_text(size = 14),
#         axis.text = element_text(size = 11))
# annotation_custom(grob)