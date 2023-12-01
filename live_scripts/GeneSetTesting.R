library(tidyverse)
library(clusterProfiler)

search_kegg_organism("mouse", by = "common_name")

shrink.d11 <- readRDS("RObjects/Shrunk_Results.d11.rds")

sigGenes <- shrink.d11 %>%
  drop_na(Entrez, padj) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  pull(Entrez)

keggRes <- enrichKEGG(gene = sigGenes, organism = 'mmu')
tab <- as_tibble(keggRes)
browseKEGG(keggRes, "mmu04612")

library(pathview)
logFC <- shrink.d11$log2FoldChange
names(logFC) <- shrink.d11$Entrez
logFC

pathview(gene.data = logFC,
         pathway.id = "mmu04612",
         species = "mmu",
         limit = list(gene=20, cpd=1))

logFC <- shrink.d11 %>%
  drop_na(padj, Entrez) %>%
  filter(padj < 0.01) %>%
  pull(log2FoldChange, Entrez)

pathview(gene.data = logFC,
         pathway.id = "mmu04659",
         species = "mmu",
         limit = list(gene=20, cpd=1))

library(org.Mm.eg.db)

sigGenes_GO <- shrink.d11 %>%
  drop_na(padj) %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
  pull(GeneID)

universe <- shrink.d11$GeneID

ego <- enrichGO(gene = sigGenes_GO,
                universe = universe,
                OrgDb = org.Mm.eg.db,
                keyType = "ENSEMBL",
                ont = "BP",
                pvalueCutoff = 0.01,
                readable = TRUE)
ego
egoTAb <- as_tibble(ego)
barplot(ego, showCategory = 20)

dotplot(ego, font.size = 14)

library(enrichplot)
ego_pt <- pairwise_termsim(ego)
emapplot(ego_pt, cex_label_category = 0.25)

library(msigdbr)

rankedGenes <- shrink.d11 %>%
  drop_na(GeneID, padj, log2FoldChange) %>%
  mutate(rank = log2FoldChange) %>%
  arrange(desc(rank)) %>%
  pull(rank, GeneID)
rankedGenes

term2gene <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, ensembl_gene)
term2name <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gs_description) %>%
  distinct()

gseaRes <- GSEA(rankedGenes,
                TERM2GENE = term2gene,
                TERM2NAME = term2name,
                pvalueCutoff = 1,
                minGSSize = 15,
                maxGSSize = 500)
gsea_tab <- as_tibble(gseaRes) %>%
  arrange(desc(abs(NES))) %>%
  top_n(10, wt=-p.adjust) %>%
  dplyr::select(-core_enrichment) %>%
  mutate(across(c("enrichmentScore", "NES"), round, digit=3)) %>%
  mutate(across(c("pvalue","p.adjust","qvalue"), scales::scientific))

gseaplot(gseaRes,
         geneSetID = "HALLMARK_INFLAMMATORY_RESPONSE",
         title = "HALLMARK_INFLAMMATORY_RESPONSE")

rankedGenes.e11 <- shrink.d11 %>%
  drop_na(GeneID, pvalue, log2FoldChange) %>%
  mutate(rank = -log10(pvalue) * sign(log2FoldChange)) %>%
  arrange(desc(rank)) %>%
  pull(rank, GeneID)

gseaRes.e11 <- GSEA(rankedGenes.e11,
                TERM2GENE = term2gene,
                TERM2NAME = term2name,
                pvalueCutoff = 1,
                minGSSize = 15,
                maxGSSize = 500)

gsea_tab.e11 <- as_tibble(gseaRes.e11) %>%
  arrange(desc(abs(NES))) %>%
  top_n(10, wt=-p.adjust) %>%
  mutate(across(c("enrichmentScore", "NES"), round, digit=3)) %>%
  mutate(across(c("pvalue","p.adjust","qvalue"), scales::scientific))

gseaplot(gseaRes.e11,
         geneSetID = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
         title = "HALLMARK_OXIDATIVE_PHOSPHORYLATION")

shrink.d33 <- readRDS("RObjects/Shrunk_Results.d33.rds")

rankedGenes.e33 <- shrink.d33 %>%
  drop_na(GeneID, pvalue, log2FoldChange) %>%
  mutate(rank = -log10(pvalue) * sign(log2FoldChange)) %>%
  arrange(desc(rank)) %>%
  pull(rank, GeneID)

gseaRes.e33 <- GSEA(rankedGenes.e33,
                    TERM2GENE = term2gene,
                    TERM2NAME = term2name,
                    pvalueCutoff = 1,
                    minGSSize = 15,
                    maxGSSize = 500)

gsea_tab.e33 <- as_tibble(gseaRes.e33) %>%
  arrange(desc(abs(NES))) %>%
  top_n(10, wt=-p.adjust) %>%
  mutate(across(c("enrichmentScore", "NES"), round, digit=3)) %>%
  mutate(across(c("pvalue","p.adjust","qvalue"), scales::scientific))




