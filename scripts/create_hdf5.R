
library(rhdf5)
library(Matrix)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(data.table)
write_hdf5 <- function(path, counts_file, organism, assembly) {
  #Read in counts file
  expr <- as.data.frame(fread(counts_file, sep = "\t", header = TRUE))
  rownames(expr) <- expr$GENE
  expr <- expr[, 2:ncol(expr)]

  #Create Files for hdf5 write
  barcodes <- colnames(expr)
  gene_name <- rownames(expr)
  matrix <- Matrix(expr, sparse = TRUE)
  feature_type <- "Gene Expression"
  #Create Ensembl Symbol group
  if (organism == "mouse") {
    gene_symbol <- mapIds(org.Mm.eg.db,
                          keys = gene_name,
                          column = "ENSEMBL",
                          keytype = "SYMBOL",
                          multiVals = "first")
    } else {
    gene_symbol <- mapIds(org.Hs.eg.db,
                          keys = gene_name,
                          column = "ENSEMBL",
                          keytype = "SYMBOL",
                          multiVals = "first")
                          }

  #Create HDF5 File
  path <- path.expand(path)
  h5createFile(path)
  group <- "matrix"
  #Write Matrix Information to HDF5 File
  h5createGroup(path, group)
  h5write(barcodes,
          file = path,
          name = paste0(group, "/barcodes"))
  h5write(matrix@x,
          file = path,
          name = paste0(group, "/data"))
  h5write(matrix@i,
          file = path,
          name = paste0(group, "/indices"))
  h5write(matrix@Dim,
          file = path,
          name = paste0(group, "/shape"))
  h5write(matrix@p,
          file = path,
          name = paste0(group, "/indptr"))

  #Write Feature info in hdf5 file
  h5createGroup(path, file.path(group, "features"))
  h5write(gene_symbol,
          file = path,
          name = paste0(group, "/features/id"))
  h5write(gene_name,
          file = path,
          name = paste0(group, "/features/name"))
  h5write(rep(feature_type, length.out = length(gene_name)),
          file = path,
          name = paste0(group, "/features/feature_type"))
  h5write("genome",
          file = path,
          name = paste0(group, "/features/_all_tag_keys"))
  h5write(rep(assembly, length.out = length(gene_name)),
          file = path,
          name = paste0(group, "/features/genome"))
}

#CREATE OPTION PARSER FOR HDF5 FUNCTION NEES
library(optparse)
option_list <- list(
  make_option(c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "count matrix .tsv format (gene by cell matrix)",
    metavar = "character"),

  make_option(c("-o", "--out"),
    type = "character",
    default = "out.hdf5",
    help = "output file to store hdf5 structure",
    metavar = "character"),

  make_option(c("-s", "--species"),
    type = "character",
    default = NULL,
    help = "genome species (mouse or human)",
    metavar = "character"),

  make_option(c("-v", "--assembly_version"),
    type = "character",
    default = NULL,
    help = "genome assembly version",
    metavar = "character")
)
#GENERATE PARSER FROM COMMDANDLINE INPUTS
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

#CALL HDF5 FUNCTION
write_hdf5(path = opt$out,
  counts_file = opt$input,
  organism = opt$species,
  assembly = opt$assembly_version)