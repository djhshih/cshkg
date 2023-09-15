# Context-specific house keeping genes

## Purpose
The purpose of this project lies in analysing the expression levels of various genes of different tissues and discovering a "context-specific" housekeeping gene. "Context-specific" refers to consistency in the gene expression across different cell types, disease and developmental stage. 

## Set up

SSH key
A public key was set and saved in /home/yerihan/.ssh/id_ed25519.pub

[Guide to set up a SSH key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)
```
ssh-keygen -t ed25519 -C "hanyeri0223@gmail.com"
```

Github
```
git clone git@github.com:hanyeri/cshkg
```
User information
```
git config --global user.email "hanyeri0223@gmail.com"
git config --global user.name "Yeri Han"
```

On Fedora 36, install ...
```
sudo dnf install rstudio-desktop
```
Developmental files-devel, R header files
```
sudo dnf install make automake gcc gcc-c++ kernel-devel
sudo dnf install R-core-devel R-java-devel libRmath-devel
```
Install R packages

For preprocessing:

* BioConductor
* BiocManager
* annotate
* tidyverse: include read_tsv() function, packages including "ragg" and "textshaping" had to be installed separately and directly on the terminal due to configuration fails.
* hgu133a.db: human genome u133a expression array annotation data
* io: read from, write to, plot to in an unified manner
```
install.packages("io")
BiocManager::install(c("annotate", "hgu133a.db"))
```

For Analysis
```
install.packages(c("ggplot2","dplyr")
```

### Remarks
Bioconductor packages could not be downloaded completely at first because R was an older 
version and the outdated packages were relocated to a backup server, taking too much time to be retrieved. Both Fedora and R are updated and biocManager packages were then be able to be downloaded. 
[Guide on updating Fedora](https://docs.fedoraproject.org/en-US/quick-docs/dnf-system-upgrade/)

[1-3-2023] loaded data in rds file to ease handling on RStudio.


## Analysis
Set R directory: 
```
[Session]-[Set Working Directory]-[To Source file location]
```
Preprocess the data
1. biogps
2. pancan
<br />[1-3-2023] explore.R: grouping samples into sample type and cancer type; removed rare samples with group frequency <10

Calculate variations within and between the groups

For each gene, between group standard deviation
```
sd= √(〖∑_g▒〖(m_g-m ̅ 〗)〗^2/(G-1))
m ̅=  1/G ∑_(g=1)^G▒m_g 
G=number of groups
m_g=mean expression of group g
N_g=number of genes

```
For a group, within group standard deviation
```
sd= √((∑_g▒〖(x_i^((g))-x ̅^((g)))〗^2 )/((∑_g▒〖N_g)〗-1))  
x ̅^((g))=mean expression in all the samples within group g
x ̅^((g))=  1/N_g  ∑_(i=1)^(N_g)▒X_i^((g)) 
N_g=number of samples
X_i^((g))=expression in sample 1 within group g

```


## Download data
run get.sh
27/01/2023 removed zip file after unzipping

run get.sh to download pancan data
14/02/2023

## Git
* add: take modified file from working directory to Git index (staging area)
* commit: creates a new revision log, you can only commit after you add
* push: add the change to the GitHub
* diff: outputs the difference between two inputs

[Guide on Git](https://rogerdudler.github.io/git-guide/)
