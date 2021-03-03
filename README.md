# scRNA-seq on microglia from wild-derived Alzheimer's mouse model

## Publication, data and shiny app
This github site stores all the scripts to reproduce the analysis in for our paper published in Cell Reports: https://doi.org/10.1016/j.celrep.2021.108739

All the raw data, processed data and analysis data (the Seurat object) have been submitted to AD Knowledge Portal (http://adknowledgeportal.synapse.org). All the contents can be found at: https://doi.org/10.7303/syn23763409

We have built a shiny app for visualizing gene expression of microglia clusters. https://wild_microglia_scrna-seq.jax.org/


## Brief Methods

Microglia were mechanically isolated from 8-9 months-old female mice from our wild-derived AD mouse panel [C57BL/6J (B6), PWK/PhJ, WSB/EiJ and CAST/EiJ strains, each carrying APP/PS1] using magnetic activated cell sorting with cd11b micro-beads and were subjected to scRNA-seq (JAX Single Cell Biology Laboratory). The whole process was performed on 4 °C. 

Data were processed using a customized pipeline (scBASE) allowing strain-specific genome alignment followed by downstream analyses using Seurat and edgeR packages. 