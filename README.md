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


## Analysis
Set R directory: 
```
[Session]-[Set Working Directory]-[To Source file location]
```
Preprocess the data


## Download data
run get.sh
27/01/2023 removed zip file after unzipping

pancan-run get.sh to download tcga data
14/02/2023


## Git
* add: take modified file from working directory to Git index (staging area)
* commit: creates a new revision log, you can only commit after you add
* push: add the change to the GitHub
* diff: outputs the difference between two inputs

[Guide on Git](https://rogerdudler.github.io/git-guide/)
