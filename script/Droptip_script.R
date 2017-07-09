## removing a certain tip
setwd("/Users/yosapol/Desktop/Last_tree/script/Final_script/")
library(ape)

source("read.newick.R") # this script is needed to read tree with singletons (like yours)

## read original tree
tr <- read.newick("new_TPS.txt")

## there are singletons, due to notational differences. Collapse these in the newick format. It won't affect the topology

tr2 <- collapse.singles(tree = tr)

### Drop this tip: AA scaffold5354 8.mRNA1|A.arabicum-Aethionema arabicum

tr3 <- drop.tip(phy = tr2, tip = c("LOC104818549|GCF_000463585.1_ASM46358v1_genomic-Cleome_hassleriana",
                                   "LOC104809000|GCF_000463585.1_ASM46358v1_genomic-Cleome_hassleriana",
                                   "LOC104818549|GCF_000463585.1_ASM46358v1_genomic-Cleome_hassleriana",
                                    "LOC104809000|GCF_000463585.1_ASM46358v1_genomic-Cleome_hassleriana",
                                   "Araha.5287s0001|Ahalleri_264_v1-Arabidopsis_halleri",
                                   "Araha.8081s0001|Ahalleri_264_v1-Arabidopsis_halleri",
                                   "Araha.20010s0001|Ahalleri_264_v1-Arabidopsis_halleri",
                                   "BjuB017318|Bju_genome-Brassica_juncea",
                                   "BniB034130|BniB_genome-Brassica_nigra",
                                   "CARUB_v10007073mg|GCF_000375325.1_Caprub1_0_genomic-Capsella_rubella",
                                   "LOC104732654|GCF_000633955.1_Cs_genomic-Camelina_sativa",
                                   "LOC10471_2471|GCF_000633955.1_Cs_genomic-Camelina_sativa",
                                   "LOC104773008|GCF_000633955.1_Cs_genomic-Camelina_sativa",
                                   "ARALYDRAFT_664809|GCF_000004255.1_v.1.0_genomic-Arabidopsis_lyrata",
                                   "AA_scaffold5354_8.mRNA1|A.arabicum-Aethionema_arabicum",
                                   "EUTSA_v10017701mg|GCF_000478725.1_Eutsalg1_0_genomic-Eutrema_salsugineum"))

             
### write this tree out as newick again

write.tree(tr3, file = "minimum_3_cluster_final.txt")

