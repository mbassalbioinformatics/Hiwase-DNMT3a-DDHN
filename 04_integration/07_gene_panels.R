#!/usr/bin/env Rscript
###############################################################################
# 07_gene_panels.R                                                   2026-10-08
#
# Gene panels for annotating DMRs, testing enrichment and choosing locus plots.
# Deliberately broad. Every gene keeps the source it came from.
#
# MSigDB (msigdbr; C2 curated, C5 GO:BP, C8 cell type, H hallmark) and curated
# lists, written out in full below with their sources:
#   HSC_stemness            Eppert 2011, Ivanova 2002, Jaatinen 2006, Park 2002, Hay 2018 (CD34+ HSC),
#                           Zheng 2018 (cord blood HSC/MPP), GO HSC homeostasis/proliferation/differentiation
#   leukaemic_stem_cell     Eppert 2011 LSC, Gal 2006, Gentles 2010, Huang AML LSC47, LSC17 (Ng 2016)
#   HOX_MEIS_self_renewal   HOX clusters A-D, MEIS1/2, PBX1/3, HOXA9-MEIS1 targets (Hess 2006, Dorsam 2004)
#   progenitors             Hay 2018 CD34+ MPP/LMPP/CLP/MEP/ERP/GRAN/MKP; Zheng 2018 CMP/GMP/MEP-like clusters
#   erythroid / megakaryocyte / granulocyte_monocyte / lymphoid
#                           Hay 2018 bone-marrow cell types, GO and Reactome differentiation sets
#   epigenetic_regulators   curated DNA methylation machinery and readers, Polycomb, COMPASS, H3K36,
#                           demethylases, acetyltransferases, cohesin, remodellers, plus GO DNA
#                           methylation/demethylation sets
#   chromatin_remodeling_GO GO:BP chromatin remodeling (broad, kept separate)
#   myeloid_drivers_CH      ELN 2022 AML, WHO/ICC myeloid, clonal haematopoiesis genes (Jaiswal 2014,
#                           Genovese 2014, Bick 2020), common myeloid NGS panels
#   splicing                curated spliceosome core and regulators, plus GO mRNA splicing via spliceosome
#   p53_targets             RNA-seq report panel, plus Hallmark p53 pathway
#   imprinted               human imprinted genes with known imprinted DMRs
#
# Output: integration/resources/gene_panels.tsv (panel, gene, source) and
#         gene_panels_provenance.tsv (panel, n genes, sources)
###############################################################################
suppressPackageStartupMessages({ library(data.table); library(msigdbr) })
source("config/paths.R")
source("03_rnaseq/00_config.R")                     # FOCUS_GENES

ms <- as.data.table(rbind(msigdbr(species = "Homo sapiens", collection = "C2"),
                          msigdbr(species = "Homo sapiens", collection = "C8"),
                          msigdbr(species = "Homo sapiens", collection = "H"),
                          msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")))
ms <- unique(ms[, .(gs_name, gene = gene_symbol)])
from_msig <- function(panel, sets) {
  miss <- setdiff(sets, ms$gs_name)
  if (length(miss)) message("  ", panel, ": not in this MSigDB release: ", paste(miss, collapse = ", "))
  ms[gs_name %in% sets, .(panel = panel, gene, source = gs_name)]
}
curated <- function(panel, genes, source) data.table(panel = panel, gene = unique(genes), source = source)

P <- list()
P$HSC_stemness <- from_msig("HSC_stemness", c(
  "EPPERT_HSC_R", "EPPERT_CE_HSC_LSC", "IVANOVA_HEMATOPOIESIS_STEM_CELL", "IVANOVA_HEMATOPOIESIS_STEM_CELL_LONG_TERM",
  "IVANOVA_HEMATOPOIESIS_STEM_CELL_SHORT_TERM", "JAATINEN_HEMATOPOIETIC_STEM_CELL_UP", "PARK_HSC_MARKERS",
  "HAY_BONE_MARROW_CD34_POS_HSC", "ZHENG_CORD_BLOOD_C6_HSC_MULTIPOTENT_PROGENITOR",
  "GOBP_HEMATOPOIETIC_STEM_CELL_HOMEOSTASIS", "GOBP_HEMATOPOIETIC_STEM_CELL_PROLIFERATION",
  "GOBP_HEMATOPOIETIC_STEM_CELL_DIFFERENTIATION"))
P$leukaemic_stem_cell <- rbind(
  from_msig("leukaemic_stem_cell", c("EPPERT_LSC_R", "GAL_LEUKEMIC_STEM_CELL_UP", "GENTLES_LEUKEMIC_STEM_CELL_UP", "HUANG_AML_LSC47")),
  curated("leukaemic_stem_cell", c("DNMT3B", "ZBTB46", "NYNRIN", "ARHGAP22", "LAPTM4B", "MMRN1", "DPYSL3", "FAM30A",
                                   "CDK6", "CPXM1", "SOCS2", "SMIM24", "EMP1", "BEX3", "CD34", "AKR1C3", "ADGRG1"),
          "LSC17 (Ng et al. 2016 Nature)"))
hox <- c(paste0("HOXA", c(1:7, 9:11, 13)), paste0("HOXB", c(1:9, 13)), paste0("HOXC", c(4:6, 8:13)), paste0("HOXD", c(1, 3, 4, 8:13)))
P$HOX_MEIS_self_renewal <- rbind(
  curated("HOX_MEIS_self_renewal", c(hox, "MEIS1", "MEIS2", "PBX1", "PBX3", "MN1", "MECOM", "HOTAIRM1"), "HOX clusters A-D, MEIS1/2, PBX1/3, MN1, MECOM"),
  from_msig("HOX_MEIS_self_renewal", c("HESS_TARGETS_OF_HOXA9_AND_MEIS1_UP", "DORSAM_HOXA9_TARGETS_UP")))
P$progenitors <- from_msig("progenitors", c(
  "HAY_BONE_MARROW_CD34_POS_MULTILIN", "HAY_BONE_MARROW_CD34_POS_LMPP", "HAY_BONE_MARROW_CD34_POS_CLP",
  "HAY_BONE_MARROW_CD34_POS_MEP", "HAY_BONE_MARROW_CD34_POS_ERP", "HAY_BONE_MARROW_CD34_POS_ERP_EARLY",
  "HAY_BONE_MARROW_CD34_POS_GRAN", "HAY_BONE_MARROW_CD34_POS_MKP", "HAY_BONE_MARROW_CD34_POS_EO_B_MAST",
  "ZHENG_CORD_BLOOD_C3_MEGAKARYOCYTE_ERYTHROID_PROGENITOR", "ZHENG_CORD_BLOOD_C9_GRANULOCYTE_MACROPHAGE_PROGENITOR",
  "ZHENG_CORD_BLOOD_C1_PUTATIVE_MEGAKARYOCYTE_PROGENITOR", "ZHENG_CORD_BLOOD_C10_MULTILYMPHOID_PROGENITOR",
  "IVANOVA_HEMATOPOIESIS_EARLY_PROGENITOR", "IVANOVA_HEMATOPOIESIS_INTERMEDIATE_PROGENITOR"))
P$erythroid <- from_msig("erythroid", c("HAY_BONE_MARROW_ERYTHROBLAST", "HAY_BONE_MARROW_EARLY_ERYTHROBLAST",
  "ZHENG_CORD_BLOOD_C4_PUTATIVE_EARLY_ERYTHROID_COMMITMENT", "GOBP_ERYTHROCYTE_DIFFERENTIATION", "HALLMARK_HEME_METABOLISM"))
P$megakaryocyte <- from_msig("megakaryocyte", c("HAY_BONE_MARROW_PLATELET", "GOBP_MEGAKARYOCYTE_DIFFERENTIATION",
  "REACTOME_FACTORS_INVOLVED_IN_MEGAKARYOCYTE_DEVELOPMENT_AND_PLATELET_PRODUCTION"))
P$granulocyte_monocyte <- from_msig("granulocyte_monocyte", c("HAY_BONE_MARROW_NEUTROPHIL", "HAY_BONE_MARROW_IMMATURE_NEUTROPHIL",
  "HAY_BONE_MARROW_MONOCYTE", "HAY_BONE_MARROW_EOSINOPHIL", "GOBP_GRANULOCYTE_DIFFERENTIATION",
  "REACTOME_TRANSCRIPTIONAL_REGULATION_OF_GRANULOPOIESIS"))
P$lymphoid <- from_msig("lymphoid", c("HAY_BONE_MARROW_CD34_POS_PRE_B", "HAY_BONE_MARROW_PRO_B", "HAY_BONE_MARROW_FOLLICULAR_B_CELL",
  "HAY_BONE_MARROW_PLASMA_CELL", "HAY_BONE_MARROW_NAIVE_T_CELL", "HAY_BONE_MARROW_CD8_T_CELL", "HAY_BONE_MARROW_NK_CELLS"))

P$epigenetic_regulators <- rbind(
  curated("epigenetic_regulators", c(
    # DNA methylation writers, erasers, cofactors and readers
    "DNMT1", "DNMT3A", "DNMT3B", "DNMT3L", "UHRF1", "UHRF2", "TET1", "TET2", "TET3", "IDH1", "IDH2", "TDG",
    "MBD1", "MBD2", "MBD3", "MBD4", "MBD5", "MBD6", "MECP2", "ZBTB33", "ZBTB4", "ZBTB38", "DPPA3", "ZFP57", "TRIM28",
    "AICDA", "APOBEC3A", "WT1", "KDM2B", "CXXC1", "FBXL19", "SALL4", "VHL",
    # Polycomb (PRC1/PRC2) and H3K27 demethylation
    "EZH1", "EZH2", "SUZ12", "EED", "JARID2", "MTF2", "PHF1", "PHF19", "AEBP2", "RBBP4", "RBBP7",
    "BMI1", "PCGF1", "PCGF2", "PCGF3", "PCGF5", "PCGF6", "RING1", "RNF2", "CBX2", "CBX4", "CBX6", "CBX7", "CBX8",
    "PHC1", "PHC2", "BCOR", "BCORL1", "KDM6A", "KDM6B", "UTY", "ASXL1", "ASXL2", "ASXL3", "BAP1",
    # COMPASS / H3K4, H3K36, H3K79, H3K9
    "KMT2A", "KMT2B", "KMT2C", "KMT2D", "KMT2E", "SETD1A", "SETD1B", "MEN1", "ASH2L", "WDR5", "RBBP5", "DPY30",
    "KDM5A", "KDM5B", "KDM5C", "KDM5D", "KDM1A", "NSD1", "NSD2", "NSD3", "SETD2", "ASH1L", "DOT1L",
    "SETDB1", "SETDB2", "EHMT1", "EHMT2", "SUV39H1", "SUV39H2", "KDM4A", "KDM4B", "KDM4C", "KDM3A", "KDM3B",
    "CBX1", "CBX3", "CBX5",
    # acetylation, cohesin, remodelling, others mutated in myeloid disease
    "CREBBP", "EP300", "KAT6A", "KAT6B", "HDAC1", "HDAC2", "SIRT1", "BRD4", "STAG1", "STAG2", "RAD21", "SMC1A",
    "SMC3", "NIPBL", "CTCF", "ARID1A", "ARID1B", "ARID2", "SMARCA2", "SMARCA4", "SMARCB1", "CHD4", "CHD8",
    "PHF6", "PRDM16", "PRDM9", "ATRX", "DAXX", "H3-3A", "H3C2", "MLLT1", "MLLT3", "AFF1", "MLLT10"),
    "curated: DNA methylation machinery and readers, Polycomb, COMPASS/H3K4, H3K36, H3K9, H3K79, acetylation, cohesin, remodellers"),
  from_msig("epigenetic_regulators", c("GOBP_DNA_METHYLATION_DEPENDENT_CONSTITUTIVE_HETEROCHROMATIN_FORMATION",
                                       "GOBP_CHROMOSOMAL_5_METHYLCYTOSINE_DNA_DEMETHYLATION_PATHWAY")))
P$chromatin_remodeling_GO <- from_msig("chromatin_remodeling_GO", "GOBP_CHROMATIN_REMODELING")

P$myeloid_drivers_CH <- curated("myeloid_drivers_CH", c(
  "ASXL1", "ASXL2", "BCOR", "BCORL1", "BRAF", "BRCC3", "CALR", "CBL", "CBLB", "CEBPA", "CHEK2", "CREBBP", "CSF3R",
  "CSNK1A1", "CUX1", "DDX3X", "DDX41", "DNMT3A", "EP300", "ETNK1", "ETV6", "EZH2", "FLT3", "GATA1", "GATA2", "GNAS",
  "GNB1", "IDH1", "IDH2", "IKZF1", "JAK2", "JAK3", "KDM6A", "KIT", "KMT2A", "KMT2C", "KMT2D", "KRAS", "MECOM", "MGA",
  "MPL", "MYC", "MYD88", "NF1", "NOTCH1", "NPM1", "NRAS", "NUP98", "PHF6", "PIGA", "PPM1D", "PRPF8", "PTPN11", "RAD21",
  "RIT1", "RUNX1", "SAMD9", "SAMD9L", "SETBP1", "SF1", "SF3A1", "SF3B1", "SH2B3", "SMC1A", "SMC3", "SRSF2", "STAG2",
  "STAT3", "STAT5B", "TET2", "TP53", "U2AF1", "U2AF2", "WT1", "ZBTB33", "ZRSR2", "ATM", "LUC7L2", "MLLT10", "PDGFRA",
  "PDGFRB", "FIP1L1", "ABL1", "ANKRD26", "SRP72", "TERT", "TERC"),
  "curated: ELN 2022 AML, WHO/ICC myeloid, clonal haematopoiesis (Jaiswal 2014, Genovese 2014, Bick 2020), myeloid NGS panels")
P$myeloid_drivers_CH <- rbind(P$myeloid_drivers_CH, curated("myeloid_drivers_CH", FOCUS_GENES$Myeloid_drivers, "RNA-seq report panel"))

spl <- intersect(c("GOBP_MRNA_SPLICING_VIA_SPLICEOSOME", "GOBP_RNA_SPLICING_VIA_TRANSESTERIFICATION_REACTIONS",
                   "GOBP_ALTERNATIVE_MRNA_SPLICING_VIA_SPLICEOSOME"), unique(ms$gs_name))
P$splicing <- rbind(
  curated("splicing", c(FOCUS_GENES$Splicing_factor, "SF1", "SF3A1", "SF3B1", "SF3B2", "SF3B3", "SF3B4", "SF3B5",
                        "SF3B6", "PHF5A", "U2AF1", "U2AF2", "SRSF1", "SRSF2", "SRSF3", "SRSF5", "SRSF6", "SRSF7",
                        "TRA2B", "HNRNPA1", "HNRNPC", "HNRNPK", "HNRNPH1", "PTBP1", "RBM10", "RBM17", "RBM25",
                        "RBM39", "PRPF3", "PRPF4", "PRPF6", "PRPF8", "PRPF31", "PRPF40B", "SNRNP200", "EFTUD2",
                        "LUC7L2", "DHX15", "DDX41", "SRRM2", "ZRSR2", "SUGP1", "QKI", "MBNL1", "CELF1", "ESRP1"),
          "curated: spliceosome core, SR and hnRNP proteins, splicing regulators mutated in myeloid disease"),
  from_msig("splicing", spl))
P$p53_targets <- rbind(curated("p53_targets", FOCUS_GENES$p53_targets, "RNA-seq report panel"),
                       from_msig("p53_targets", "HALLMARK_P53_PATHWAY"))
P$imprinted <- curated("imprinted", c(
  "PEG3", "ZIM2", "PEG10", "SGCE", "MEST", "MESTIT1", "H19", "IGF2", "KCNQ1", "KCNQ1OT1", "CDKN1C", "PHLDA2",
  "SLC22A18", "SNRPN", "SNURF", "MKRN3", "MAGEL2", "NDN", "NPAP1", "SNHG14", "UBE3A", "ATP10A", "NNAT", "BLCAP",
  "GNAS", "GNAS-AS1", "GRB10", "DLK1", "MEG3", "RTL1", "DIO3", "PLAGL1", "HYMAI", "ZNF331", "ZDBF2", "L3MBTL1",
  "NAP1L5", "PEG13", "KCNK9", "DIRAS3", "INPP5F", "FAM50B", "ZNF597", "NAA60", "MIMT1", "ZNF264", "TP73", "GPR1",
  "WT1", "DLGAP2", "AIM1", "LIN28B", "RB1", "HM13", "ERAP2"),
  "curated: human imprinted genes with known imprinted DMRs (geneimprint.com, Monk et al. 2014, Court et al. 2014)")

pan <- unique(rbindlist(P), by = c("panel", "gene", "source"))
fwrite(pan, file.path(RES, "gene_panels.tsv"), sep = "\t")
prov <- pan[, .(genes = uniqueN(gene), sources = paste(unique(source), collapse = "; ")), by = panel]
fwrite(prov, file.path(RES, "gene_panels_provenance.tsv"), sep = "\t")
print(prov[, .(panel, genes, n_sources = lengths(strsplit(sources, "; ")))])
