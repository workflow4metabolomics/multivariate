#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = F)
script.path <- sub("--file=","",args[grep("--file=",args)])
library(RUnit)
library(ropls)

#############
# CONSTANTS #
#############

exaDirOutC <- file.path(dirname(script.path), 'output')

###################
# LOAD DATA FRAME #
###################

load.df <- function(file) {
	file.exists(file) || stop(paste0("No output file \"", file ,"\"."))
	return(read.table(file = file, header = TRUE, sep = "\t", stringsAsFactors = FALSE))
}

#############################
# CALL MULTIVARIATE WRAPPER #
#############################

call.multivariate.wrapper <- function(params) {

	# Set program path
	prog <- file.path(dirname(script.path), '..', 'multivariate_wrapper.R')

	# Set arguments
	args <- NULL
	for (a in names(params))
		args <- c(args, a, params[[a]])

	# Call
	call <- paste(c(prog, args), collapse = ' ')
	retcode <- system(call)
	if (retcode != 0)
		stop("Error when running multivariate_wrapper.R.")

	# Get output
	out <- list()
	if ('ropls_out' %in% names(params)) {
		load(file = params[['ropls_out']])
		out[['ropLs']] <- ropLs
	}
	if ('sampleMetadata_out' %in% names(params))
		out[['samDF']] <- load.df(params[['sampleMetadata_out']])
	if ('variableMetadata_out' %in% names(params))
		out[['varDF']] <- load.df(params[['variableMetadata_out']])

	return(out)
}

########
# MAIN #
########

## Create output folder
##---------------------
file.exists(exaDirOutC) || dir.create(exaDirOutC)

# Define data sets
tesArgLs <- list(input_pca = c(respC = "none",
                     predI = "NA",
                     orthoI = "0",
                     testL = "FALSE"),
                 input_pcaGender = c(respC = "none",
                     predI = "NA",
                     orthoI = "0",
                     testL = "FALSE",
                     parMahalC = "gender"),
                 input_plsdaGender = c(respC = "gender",
                     predI = "NA",
                     orthoI = "0",
                     testL = "FALSE"),
                 input_oplsAge = c(respC = "age",
                     predI = "1",
                     orthoI = "1",
                     testL = "FALSE"),
                 input_oplsdaGender = c(respC = "gender",
                     predI = "1",
                     orthoI = "1",
                     testL = "FALSE"),
                 sacurine_pca = c(respC = "none",
                     predI = "NA",
                     orthoI = "0",
                     testL = "FALSE"),
                 sacurine_pcaGender = c(respC = "none",
                     predI = "NA",
                     orthoI = "0",
                     testL = "FALSE",
                     parMahalC = "gender"),
                 sacurine_plsAge = c(respC = "age",
                     predI = "NA",
                     orthoI = "0",
                     testL = "FALSE"),
                 sacurine_plsdaGender = c(respC = "gender",
                     predI = "NA",
                     orthoI = "0",
                     testL = "FALSE",
                     .chkC = "checkEqualsNumeric(getSummaryDF(outLs[['ropLs']])[, 'Q2(cum)'], 0.584, tolerance = 1e-3)"),
                 sacurineTest_pls = c(respC = "age",
                     predI = "2",
                     orthoI = "0",
                     testL = "TRUE",
                     .chkC = "checkEqualsNumeric(outLs[['samDF']][181, 'age_PLS_predictions'], 40.82252, tolerance = 1e-5)"),
                 sacurineTest_opls = c(respC = "age",
                     predI = "1",
                     orthoI = "2",
                     testL = "TRUE",
                     .chkC = "checkEqualsNumeric(outLs[['samDF']][181, 'age_OPLS_predictions'], 40.28963, tolerance = 1e-5)"),
                 sacurineTest_plsda = c(respC = "gender",
                     predI = "2",
                     orthoI = "0",
                     testL = "TRUE",
                     .chkC = "checkEquals(outLs[['samDF']][181, 'gender_PLSDA_predictions'], 'F')"),
                 sacurineTest_oplsda = c(respC = "gender",
                     predI = "1",
                     orthoI = "1",
                     testL = "TRUE",
                     .chkC = "checkEquals(outLs[['samDF']][181, 'gender_OPLSDA_predictions'], 'F')"),
                 sacurine_oplsAge = c(respC = "age",
                     predI = "1",
                     orthoI = "1",
                     testL = "FALSE",
                     .chkC = "checkEqualsNumeric(outLs[['varDF']][1, 'age_OPLS_VIP_ortho'], 0.3514378, tolerance = 1e-7)"),
                 sacurine_oplsdaGender = c(respC = "gender",
                     predI = "1",
                     orthoI = "1",
                     testL = "FALSE"),
                 example1_plsda = c(respC = "traitment",
                     predI = "3",
                     orthoI = "0",
                     testL = "FALSE",
                     .chkC = "checkEqualsNumeric(nrow(outLs[['ropLs']]@modelDF), 3, tolerance = 0.5)"),
                 example2_pca = c(respC = "none",
                     predI = "NA",
                     orthoI = "0",
                     testL = "FALSE"))

# Add file information for each dataset
for(tesC in names(tesArgLs))
    tesArgLs[[tesC]] <- c(tesArgLs[[tesC]],
                          dataMatrix_in = file.path(unlist(strsplit(tesC, "_"))[1], "dataMatrix.tsv"),
                          sampleMetadata_in = file.path(unlist(strsplit(tesC, "_"))[1], "sampleMetadata.tsv"),
                          variableMetadata_in = file.path(unlist(strsplit(tesC, "_"))[1], "variableMetadata.tsv"),
                          sampleMetadata_out = file.path(exaDirOutC, "sampleMetadata.tsv"),
                          variableMetadata_out = file.path(exaDirOutC, "variableMetadata.tsv"),
                          ropls_out = file.path(exaDirOutC, "ropls.bin"),
                          figure = file.path(exaDirOutC, "figure.pdf"),
                          information = file.path(exaDirOutC, "information.txt"))

# Run tests on each dataset
for(tesC in names(tesArgLs)) {
    print(tesC)
	args <- as.list(tesArgLs[[tesC]])
	args[['.chkC']] <- NULL
    outLs <- call.multivariate.wrapper(args)
    if(".chkC" %in% names(tesArgLs[[tesC]]))
        stopifnot(eval(parse(text = tesArgLs[[tesC]][[".chkC"]])))
}

message("Checks successfully completed")
