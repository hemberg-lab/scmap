## scmap - A tool for unsupervised projection of single cell RNA-seq data

Single-cell RNA-seq (scRNA-seq) is widely used to investigate the composition of complex tissues since the technology allows researchers to define cell-types using unsupervised clustering of the transcriptome. However, due to differences in experimental methods and computational analyses, it is often challenging to directly compare the cells identified in two different experiments. Here, we present scmap, a method (source available at https://github.com/hemberg-lab/scmap and the application can be run from http://www.hemberg-lab.cloud/scmap) for projecting cells from a scRNA-seq experiment on to the cell-types identified in a different experiment.

### Questions

__Q__: How to install/run __scmap__?  
__A__: In your R session:
```
install.packages("devtools")
devtools::install_github("hemberg-lab/scmap")
```

__Q__: Where can I report bugs, comments, issues or suggestions?  
__A__: Please use [this page](https://github.com/hemberg-lab/scmap/issues).

__Q__: Is __scmap__ published?  
__A__: Not yet, but a copy of __scmap__ manuscript is available on [bioRxiv](http://doi.org/10.1101/150292).

__Q__: What is __scmap__ licence?  
__A__: GPL-3
