# Download Touchstone Level 2 GEX delta dataset
wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5FLevel2%5FGEX%5Fdelta%5Fn49216x978%2Egctx%2Egz'
gunzip GSE92742_Broad_LINCS_Level2_GEX_delta_n49216x978.gctx.gz

# Download Touchstone information dataset
## sample information
wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Finst%5Finfo%2Etxt%2Egz'
gunzip GSE92742_Broad_LINCS_inst_info.txt.gz

## gene and landmark information
wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fgene%5Finfo%5Fdelta%5Flandmark%2Etxt%2Egz'
gunzip GSE92742_Broad_LINCS_gene_info_delta_landmark.txt.gz
