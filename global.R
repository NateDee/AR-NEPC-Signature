library(data.table)
library(plyr)
if(Sys.getenv('SHINY_PORT') == "") {
options(shiny.maxRequestSize=1000*1024^2)
}
dbgap_fpkm = as.data.frame(fread("./datasets/DBGAP_AR_NEPC_Genes_log2fpkm_matrix.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
dbgap_scores= as.data.frame(fread("./datasets/DBGAP_AR_NEPC_Scores_matrix.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE))
nepc_ref = read.table("./datasets/NEPC_reference_vector.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
ar_ref = read.table("./datasets/AR_reference_vector.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
