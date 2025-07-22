# Context-specific house keeping genes

## Purpose

The purpose of this project lies in analysing the expression levels of various genes of different tissues and discovering a "context-specific" housekeeping gene. "Context-specific" refers to consistency in the gene expression across different cell types, disease and developmental stage. 

## Set up

### SSH

SSH key
A public key can be set and saved in /<user_path>/.ssh/id_ed25519.pub

[Guide to set up a SSH key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)
```
ssh-keygen -t ed25519 -C "<user_email_address>"
```

Github
```
git clone git@github.com:yerimeeei/cshkg
```
User information
```
git config --global user.email "<user_email_address>"
git config --global user.name "<user_name>"
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
4. touchstone

## Download and preprocess data

`cd` to the `data` subdirectory.

Within a dataset directory, run `./get.sh`.

Prepare the data by running `Rscript preprocess.R`.

## Analysis

Set R working directory to be the same as the script.
```
[Session] -> [Set Working Directory] -> [To Source file location]
```

Then, run `analyze.R` in specific dataset directory.

## Mathematcis

We calculate variations within and between the groups as follows.

For each gene, between-group standard deviation `\sigma_{bw}` is calculated as
```
$`\sigma_{bw} = \sqrt{ (\sum_k^K m_k - \bar{m})^2 / (K - 1) }`$
$`\bar{m} = \frac{1}{K} \sum_k^K m_k`$
```
where `K` is the number of groups,
`m_k` is the mean expression within group `k`,
and `m` is the overall mean expression of samples.

The within group standard deviation `\sigma_{wi}` is calculated by
```
$`\sigma_{wi} = \sqrt{ (\sum_k^K (\sum_i^N_k (x_i^{(k)} - \bar{x}^{(k)})^2)) / ( \sum_k^K N_k - 1 )  }`$
$`\bar{x}^{(k)} = \frac{1}{N_k} \sum_i^{N_k} x_i^{(k)}`$
```
where `N_k` is the number of samples within group `k` where k ∈ K,
and `x_i^{(k)}` is the expression of sample `i` within group `k`.

## Git

* diff: outputs the difference between two inputs
* status: displays the state of the working directory and the staging area
* add: take modified file from working directory to Git index (staging area)
* commit: creates a new revision log, you can only commit after you add
* push: add the change to the GitHub
* pull: updates the local repository to include changes already included in the remote repository

[Guide on Git](https://rogerdudler.github.io/git-guide/)

## Troubleshooting

### Bioconductor

Bioconductor packages could not be downloaded completely at first because R was an older 
version and the outdated packages were relocated to a backup server, taking too much time to be retrieved.  
Both Fedora and R are updated and biocManager packages were then be able to be downloaded. 
[Guide on updating Fedora](https://docs.fedoraproject.org/en-US/quick-docs/dnf-system-upgrade/)

