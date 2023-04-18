###########################################################################
## Get functions for checkpoint 2: identification of malignant subclones ##
############################################################################
#
# Written by Kristen Michelle Nader <kristen.nader@helsinki.fi>
# and Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi> , April 2021
#
# Functions on this page:
# get_gene_annotation,get_normal_cells,run_infercnv
#

lapply(c("dplyr","Seurat","HGNChelper","biomaRt","infercnv2"), library, character.only = T)
