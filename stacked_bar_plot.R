# Plot
library(data.table)
library(readr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library("RColorBrewer")


write.dataset <- function(dataset, path, sep='\t', row.names=F){
  write.table(x = dataset, file = path, sep = sep, row.names = row.names, quote = F)
}

NoLegend <- function (...){
  no.legend.theme <- theme(legend.position = "none", validate = TRUE,
                           ...)
  return(no.legend.theme)
}

NoXaxis <- function(...){
  no.x.axis.theme <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                           axis.ticks.x=element_blank(), ...)
  return(no.x.axis.theme)
}

AdjThemeY <- function(text_size=10, label_size=13, text_angle=0, label_angle=90, ...){
  adj.theme <- theme(plot.margin = unit(c(.05,.05,.05,.05), "cm"),
                     axis.text.y = element_text(angle = text_angle, size=text_size),
                     axis.title.y = element_text(angle = label_angle, size=label_size),
                     panel.grid.major.y = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank())
  
  return(adj.theme)
}

exp_pairs <- c('BC2_BC3', 'BC4_BC5')
KIR.subset$bcexp_pair_or_sin <- ifelse(!grepl('_', KIR.subset$barcode) | KIR.subset$barcode %in% exp_pairs,
                                       T, F)
KIR.subset$bcexp_pair_or_sin[is.na(KIR.subset$barcode)] <- NA

KIR.subset$barcode[(KIR.subset$HTO_classification == "S9") & 
                     (KIR.subset$barcode=="BC14_BC2") &
                     (!is.na(KIR.subset$barcode))] <-"BC14"

KIR.subset$barcode[(KIR.subset$HTO_classification == "S9") & 
                     (KIR.subset$barcode=="BC14_BC7") &
                     (!is.na(KIR.subset$barcode))] <-"BC14"

KIR.subset$barcode[(KIR.subset$HTO_classification == "S9") & 
                     (KIR.subset$barcode=="BC14_BC6") &
                     (!is.na(KIR.subset$barcode))] <-"BC14"

KIR.subset$barcode[(KIR.subset$HTO_classification == "S9") & 
                     (KIR.subset$barcode=="BC14_BC5") &
                     (!is.na(KIR.subset$barcode))] <-"BC14"


## Buffy coat 357

condition_clone_bar_BC357 <- KIR.subset@meta.data %>% 
  dplyr::select(c(topb_clones_BC357, condition)) %>% 
  filter(!is.na(topb_clones_BC357)) %>% 
  mutate(condition = factor(condition)) %>% 
  group_by(condition, topb_clones_BC357) %>% 
  summarise(clone_freq=n()) %>% 
  ggplot(aes(x=topb_clones_BC357, y=clone_freq, fill=condition)) +
  geom_bar(stat="identity") +
  labs(x="beta clone ID",y="",size='GEMs',color='GEMs') +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  labs(y="Condition") 


barcode_clone_bar_BC357 <- KIR.subset@meta.data %>% 
  filter(bcexp_pair_or_sin) %>%
  dplyr::select(c(topb_clones_BC357, barcode)) %>% 
  filter(!is.na(topb_clones_BC357)) %>% 
  mutate(barcode = factor(barcode)) %>% 
  group_by(barcode, topb_clones_BC357) %>% 
  summarise(clone_freq=n()) %>% 
  ggplot(aes(x=topb_clones_BC357, y=clone_freq, fill=barcode)) +
  scale_fill_manual(values=brewer.pal(n = 12, name = "Paired")) +
  geom_bar(stat="identity") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  labs(y="Barcode") 


hash_clone_bar_BC357 <- KIR.subset@meta.data %>% 
  dplyr::select(c(topb_clones_BC357, HTO_classification)) %>% 
  filter(!is.na(topb_clones_BC357)) %>% 
  mutate(HTO_classification = factor(HTO_classification)) %>% 
  group_by(HTO_classification, topb_clones_BC357) %>% 
  summarise(clone_freq=n()) %>% 
  ggplot(aes(x=topb_clones_BC357, y=clone_freq, fill=HTO_classification)) +
  scale_fill_manual(values=brewer.pal(n = 5, name = "Set2"))+
  geom_bar(stat="identity") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  labs(y="Hashing") 

OverviewPlot_BC357 <- cowplot::plot_grid(condition_clone_bar_BC357 + NoLegend() + NoXaxis() + AdjThemeY(label_size=10),
                                   barcode_clone_bar_BC357 + NoLegend() + NoXaxis() + AdjThemeY(label_size=10),
                                   hash_clone_bar_BC357 + NoLegend() + AdjThemeY(label_size=10),
                                   align = "v", ncol = 1, rel_heights = c(.2, .2, .3), rel_widths = c(1, 1))
OverviewPlot_BC357


legend_ccb_BC357 <- cowplot::get_legend(condition_clone_bar_BC357)
legend_bcb_BC357 <- cowplot::get_legend(barcode_clone_bar_BC357 )
legend_hcb_BC357 <- cowplot::get_legend(hash_clone_bar_BC357 )

legends_mutbar_BC357 <- cowplot::plot_grid(legend_ccb_BC357, legend_bcb_BC357, legend_hcb_BC357,
                                     align = "v", ncol = 3, rel_widths = c(1, 1))

legends_mutbar_BC357


OverviewPlot_BC357 <- cowplot::plot_grid(OverviewPlot_BC357, legends_mutbar_BC357, align = "v", ncol = 1, rel_heights = c(.85, .15), rel_widths = c(1, 1))
OverviewPlot_BC357





###################################################

## Buffy coat 351

condition_clone_bar_BC351 <- KIR.subset@meta.data %>% 
  dplyr::select(c(topb_clones_BC351, condition2)) %>% 
  filter(!is.na(topb_clones_BC351)) %>% 
  mutate(condition2 = factor(condition2)) %>% 
  group_by(condition2, topb_clones_BC351) %>% 
  summarise(clone_freq=n()) %>% 
  ggplot(aes(x=topb_clones_BC351, y=clone_freq, fill=condition2)) +
  geom_bar(stat="identity") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  labs(y="Condition") 


barcode_clone_bar_BC351 <- KIR.subset@meta.data %>% 
  filter(bcexp_pair_or_sin) %>%
  dplyr::select(c(topb_clones_BC351, barcode)) %>% 
  filter(!is.na(topb_clones_BC351)) %>% 
  mutate(barcode = factor(barcode)) %>% 
  group_by(barcode, topb_clones_BC351) %>% 
  summarise(clone_freq=n()) %>% 
  ggplot(aes(x=topb_clones_BC351, y=clone_freq, fill=barcode)) +
  scale_fill_manual(values=brewer.pal(n = 12, name = "Paired")) +
  geom_bar(stat="identity") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  labs(y="Barcode") 


hash_clone_bar_BC351 <- KIR.subset@meta.data %>% 
  dplyr::select(c(topb_clones_BC351, HTO_classification)) %>% 
  filter(!is.na(topb_clones_BC351)) %>% 
  mutate(HTO_classification = factor(HTO_classification)) %>% 
  group_by(HTO_classification, topb_clones_BC351) %>% 
  summarise(clone_freq=n()) %>% 
  ggplot(aes(x=topb_clones_BC351, y=clone_freq, fill=HTO_classification)) +
  scale_fill_manual(values=brewer.pal(n = 5, name = "Set2"))+
  geom_bar(stat="identity") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  labs(y="Hashing") 

OverviewPlot_BC351 <- cowplot::plot_grid(condition_clone_bar_BC351 + NoLegend() + NoXaxis() + AdjThemeY(label_size=10),
                                   barcode_clone_bar_BC351 + NoLegend() + NoXaxis() + AdjThemeY(label_size=10),
                                   hash_clone_bar_BC351 + NoLegend() + AdjThemeY(label_size=10),
                                   align = "v", ncol = 1, rel_heights = c(.2, .2, .3), rel_widths = c(1, 1))
OverviewPlot_BC351


legend_ccb_BC351 <- cowplot::get_legend(condition_clone_bar_BC351)
legend_bcb_BC351 <- cowplot::get_legend(barcode_clone_bar_BC351)
legend_hcb_BC351 <- cowplot::get_legend(hash_clone_bar_BC351)

legends_mutbar_BC351 <- cowplot::plot_grid(legend_ccb_BC351, legend_bcb_BC351, legend_hcb_BC351,
                                           align = "v", ncol = 3, rel_widths = c(1, 1))
legends_mutbar_BC351


OverviewPlot_BC351 <- cowplot::plot_grid(OverviewPlot_BC351, legends_mutbar_BC351, align = "v", ncol = 1, rel_heights = c(.85, .15), rel_widths = c(1, 1))
OverviewPlot_BC351









