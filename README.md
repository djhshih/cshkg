# Context-specific house keeping genes

## Purpose

The purpose of this project lies in analysing the expression levels of various genes of different tissues and discovering a "context-specific" housekeeping gene. "Context-specific" refers to consistency in the gene expression across different cell types, disease and developmental stage. 

## Set up

### SSH

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

### Rstudio

On Fedora 36, install ...
```
sudo dnf install rstudio-desktop
```
Developmental files-devel, R header files
```
sudo dnf install make automake gcc gcc-c++ kernel-devel
sudo dnf install R-core-devel R-java-devel libRmath-devel
```

### R packages

* io: read from, write to, plot to in an unified manner
* ggplot2
* dplyr
* annotate
* hgu133a.db: human genome u133a expression array annotation data

Bioconductor packages
```
install.packages("BiocManager")
BiocManager::install(c("annotate", "hgu133a.db"))
```

CRAN packages
```
install.packages(c("io", "ggplot2", "dplyr"))
```

## Datasets

This repo is organized by different datasets:

1. biogps
2. pancan
3. gtex

## Download and preprocess data

`cd` to the `data` subdirectory.

Within a dataset directory, run `./get.sh`.

Prepare the data by running `Rscript preprocess.R`.

## Analysis

Set R working directory to be the same as the script.
```
[Session] -> [Set Working Directory] -> [To Source file location]
```

Then, run `analyze.R` in pancan.

## Mathematcis

We calculate variations within and between the groups as follows.

For each gene, between-group standard deviation `\sigma_{bw}` is calculated as
```
\sigma_{bw} = \sqrt{ (\sum_g^G m_g - \bar{m})^2 / (G - 1) }
\bar{m} = \frac{1}{G} \sum_g^G m_g
```
where `G` is the number of groups,
and `m_g` is the mean expression within group `g`.

The within group standard deviation `\sigma_{wi}` is calculated by
```
\sigma_{wi} = \sqrt{ (\sum_g^G ((x_i^{(g)} - \bar{x}^{(g)})^2) / ( \sum_g^G N_g - 1 )  }
\bar{x}^{(g)} = \frac{1}{N_g} \sum_i^{N_g} x_i^{(g)}
```
where `N_g` is the number of samples within group `g`, and
`x_i^{(g)}` is the expression of sample `i` within group `g`.

## Git

* add: take modified file from working directory to Git index (staging area)
* commit: creates a new revision log, you can only commit after you add
* push: add the change to the GitHub
* diff: outputs the difference between two inputs

[Guide on Git](https://rogerdudler.github.io/git-guide/)

## Troubleshooting

### Bioconductor

Bioconductor packages could not be downloaded completely at first because R was an older 
version and the outdated packages were relocated to a backup server, taking too much time to be retrieved.  
Both Fedora and R are updated and biocManager packages were then be able to be downloaded. 
[Guide on updating Fedora](https://docs.fedoraproject.org/en-US/quick-docs/dnf-system-upgrade/)

