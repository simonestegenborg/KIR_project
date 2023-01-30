library(alakazam)

bf1 <- kir@meta.data %>% 
  filter(buffy_coat == 'BC351') %>%
  filter(!is.na(TRB_clone)) %>%
  mutate(groups = recode(HTO_classification, 'S1'='S1', 'S2'='S2_S3', 'S3'='S2_S3', 'S5'='S5'))

curve1 <- estimateAbundance(bf1, group="groups", ci=0.95, nboot=100, clone="TRB_clone", )

bf2 <- kir@meta.data %>% 
  filter(buffy_coat == 'BC357') %>%
  filter(!is.na(TRB_clone)) %>%
  mutate(groups = recode(HTO_classification, 'S6'='S6_S7', 'S7'='S6_S7', 'S8'='S8', 'S10'='S10'))

curve2 <- estimateAbundance(bf2, group="groups", ci=0.95, nboot=100, clone="TRB_clone", )

plot(curve1, legend_title="Sample") + plot(curve2, legend_title="Sample")
