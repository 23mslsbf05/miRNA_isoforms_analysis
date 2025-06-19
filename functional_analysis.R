library(clusterProfiler)
library(ReactomePA)
library(org.Rn.eg.db)
library(enrichplot)
library(ggplot2)

# Example gene symbols
gene_symbols <- c("Akr1c14", "Anxa2", "Bgn", "C4a", "C7h8orf33", "Cdh5", "Cpeb4", "Cplane2", "Cplx2", "Cpm", "Cpn1", "Cpne5", "Cpsf2", "Cptp", "Cxcl16", "Cytl1", "Dcdc2c", "Ddr2", "Des", "Drd1", "Ecel1", "Fgf14", "Fmod", "Gjb2", "Igf2r", "Islr", "Loxhd1", "Lrrc32", "Lum", "Myl9", "Nexn", "Nid1", "Ptgfr", "Ptpn14", "Ret", "Retreg3", "Slc16a12", "Slc22a6", "Slc22a8", "Synpo2l", "Tbx18", "Tgm2", "Thbs2", "Wnt6")

# Convert to Entrez IDs
entrez_ids <- bitr(gene_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Rn.eg.db")$ENTREZID

# GO Enrichment - BP
go_enrich_BP <- enrichGO(gene = entrez_ids, OrgDb = org.Rn.eg.db, keyType = "ENTREZID", 
                         ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.05, readable = TRUE)

# View results
head(go_enrich_BP)

# Save plot
#p1 <- dotplot(go_enrich_BP, showCategory = 20) + ggtitle("GO Biological Process Enrichment")
#ggsave("GO_BP_enrichment.png", p1, width = 10, height = 8, dpi = 300)

# GO Enrichment - MF
go_enrich_MF <- enrichGO(gene = entrez_ids, OrgDb = org.Rn.eg.db, keyType = "ENTREZID", 
                         ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.05, readable = TRUE)
p2 <- dotplot(go_enrich_MF, showCategory = 20) + ggtitle("GO Molecular Function Enrichment")
ggsave("GO_MF_enrichment.png", p2, width = 10, height = 8, dpi = 300)

# View results
head(go_enrich_MF)

# GO Enrichment - CC
go_enrich_CC <- enrichGO(gene = entrez_ids, OrgDb = org.Rn.eg.db, keyType = "ENTREZID", 
                         ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.05, readable = TRUE)
p3 <- dotplot(go_enrich_CC, showCategory = 20) + ggtitle("GO Cellular Component Enrichment")
ggsave("GO_CC_enrichment.png", p3, width = 10, height = 8, dpi = 300)

# View results
head(go_enrich_CC)

# KEGG Enrichment
kegg_enrich <- enrichKEGG(gene = entrez_ids, organism = 'rno', pvalueCutoff = 0.05)
kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
p4 <- dotplot(kegg_enrich, showCategory = 20) + ggtitle("KEGG Pathway Enrichment")
ggsave("KEGG_enrichment.png", p4, width = 10, height = 8, dpi = 300)

# View results
head(kegg_enrich)

# Reactome Enrichment
reactome_enrich <- enrichPathway(gene = entrez_ids, organism = "rat", 
                                 pvalueCutoff = 0.05, readable = TRUE)
p5 <- dotplot(reactome_enrich, showCategory = 20) + ggtitle("Reactome Pathway Enrichment")
ggsave("Reactome_enrichment.png", p5, width = 10, height = 8, dpi = 300)

# View results
head(reactome_enrich)
