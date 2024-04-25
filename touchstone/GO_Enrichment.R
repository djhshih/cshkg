# Load required libraries for the GO enrichment analysis.
library(io)
library(dplyr)
library(gprofiler2)

load("out/candidate_cshkg.rds")
load("out/common_hkg.rds")

# Enrichment analysis
# 1. Candidate cshkg
gp <- gost(query = list("candidate_cshkg" = hkg$gene), organism = "hsapiens", 
           significant = FALSE)

# Order genes by mean expression
hkg_ordered <- hkg[order(hkg$mean, decreasing = TRUE), ]
gp_ordered <- gost(query = list("candidate_cshkg" = hkg_ordered$gene), organism = "hsapiens", ordered_query = TRUE,
                   significant = FALSE)

# Visualization of results
p <- gostplot(gp, interactive = FALSE)
p <- gostplot(gp_ordered, interactive = FALSE)

# highlight 10 GO terms important for cellular proliferation
# qdraw(
#   publish_gostplot(p, highlight_terms = c("GO:0005737", "GO:0043229", 
#                                           "GO:0005515", "GO:0016020", "GO:0051171", 
#                                           "GO:0031982", "GO:0003723", "GO:0097708", 
#                                           "GO:0031410", "GO:0008283", "GO:0048306")),
#   width = 10, height = 9,
#   file = "plots/GO_cshkg.png"
# )
qdraw(
  publish_gostplot(p, highlight_terms = c("GO:0001725", "GO:0097517",
                                          "GO:0042641", "GO:0030906", "GO:0032432",
                                          "GO:0097422", "GO:0019842", "GO:0008483")),
  width = 10, height = 9,
  file = "plots/GO_cshkg.png"
)

# 2. known common hkg
# select significant terms
gp2 <- gost(query = list("known common hkgs" = common_hkg$gene),
                 organism = "hsapiens", significant = TRUE)
p2 <- gostplot(gp2, interactive = FALSE)
qdraw(
  publish_gostplot(p2, highlight_terms = c(gp2$result$term_id)),
  width = 15, height = 9,
  file = "plots/GO_known_significant.png"
)

# 3. candidate cshkgs vs known common hkgs
# Analyze multiple gene lists
gp_multi <- gost(query = list("candidate cshkgs" = hkg$gene,
                              "known common hkgs" = common_hkg$gene),
                 organism = "hsapiens", significant = FALSE)
p_multi <- gostplot(gp_multi, interactive = FALSE)
qdraw(
  publish_gostplot(p_multi, highlight_terms = c("GO:0005737", "GO:0043229", "GO:0005515",
                                                "GO:0016020", "GO:0051171", "GO:0031982",
                                                "GO:0003723", "GO:0097708", "GO:0031410",
                                                "GO:0008283", "GO:0046099", "GO:0052657",
                                                "REAC:R-HSA-9735804")),
  width = 15, height = 9,
  file = "plots/GO_multi.png"
)
