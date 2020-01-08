#!/usr/bin/env Rscript

Install_Multiples_Packages <- function(packages) {
    pack <- packages[!(packages %in% installed.packages()[,'Package'])];
    if (length(pack)) {
        install.packages(pack, repos = 'https://cran.rstudio.com/')
    }
    
    for (package_i in packages) {
        suppressPackageStartupMessages(library(package_i, character.only = TRUE, quietly = TRUE))
    }
    
}

Install_Multiples_Packages(c('optparse', 'parallel', "tools", "data.table", "tidyr", "dplyr"))

#suppressPackageStartupMessages(library("tools"))
#suppressPackageStartupMessages(library("parallel"))
#suppressPackageStartupMessages(library("optparse"))
#suppressPackageStartupMessages(library("data.table"))
#suppressPackageStartupMessages(library("tidyr"))

#where [GB] is an upper bound on the memory pool in gigabytes (e.g. –Xmx50g), and [arguments] is a space separated list of parameter values, each having the format parameter=value.
option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = "phased_haplotype.vcf",
                help = "The vcf file containing phased haplotypes [default %default]",
                dest = "inputFile"),
    make_option(c("-f", "--phenotype"), type = "character", default = "phenotype.txt",
                help = "phenotype file [default %default]",
                dest = "phenotypeFile"),
    make_option(c("-o", "--output"), type = "character", default = 'outputFile.txt',
                help = "Name of the result file [default %default]",
                dest = "output"),
    make_option(c("-r", "--inputFolder"), type = "character", default = "00-GenotypedFiles",
                help = "Directory where the sequence data is stored [default %default]",
                dest = "inputFolder"),
    make_option(c("-a", "--map"), type = "character", default = "genetic_map.ped",
                help = "Specifies a PLINK format genetic map with cM units [default %default]",
                dest = "map"),
    make_option(c("-b", "--haplotypeFolder"), type = "character", default = '01-haplotypeFolder',
                help = "Directory where to store the mapping results [default %default]",
                dest = "haplotypeFolder"),
    make_option(c('-w', '--window'), type  =  'character', default  = '5:6567656-6675656',
                help = 'chrom:start-end specifies a chromosome interval: chrom is the CHROM field in the input VCF file and start and end] are the starting and ending positions. The entire chromosome, the beginning, or the end may be specified by chrom; chrom:-end, and chrom:start- respectively [default %default]',
                dest = 'window'),
    make_option(c("-e", "--excludeSamples"), type = "character", default = "samplesToExclude.txt",
                help = "Specifies a file containing samples (one sample identifier per line) to be excluded from the analysis. [default %default]",
                dest = "excludeSamples"),
    make_option(c("-E", "--excludeMarkers"), type = "character", default = "markesToExclude.txt",
                help = "Specifies a file containing markers (one marker per line) to be excluded from the analysis. Each line of the file can be either an identifier from a VCF record’s ID field or a genomic coordinate in the format: CHROM:POS. [default %default]",
                dest = "excludeMarkers"),
    make_option(c("-p", "--processors"), type = "integer", default = 8,
                help = "number of processors to use [defaults %default]",
                dest = "procs"),
    make_option(c("-m", "--memmory"), type = "integer", default = 50,
                help = "here [GB] is an upper bound on the memory pool in gigabytes (e.g. -m 50g). [defaults %default]",
                dest = "memmory"),
    make_option(c("-q", "--sampleprocs"), type = "integer", default = 2,
                help = "number of samples to process at time [default %default]",
                dest = "mprocs"),
    make_option(c("-x", "--external"), action  =  'store', type  =  "character", default = 'FALSE',
                help = "A space delimeted file with a single line contain several external parameters from STAR [default %default]",
                dest = "externalParameters")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list = option_list, usage = paste('%prog', '-i', 'input_file.vcf.gz', '-o', 'output_file.dat', '-f', 'phenotypic_file.txt', '[options]')))



if (!file.exists('beagle')) {
    dir.create(file.path('beagle'), recursive = TRUE, showWarnings = FALSE)
    system(paste('wget https://faculty.washington.edu/browning/beagle/beagle.21Sep19.ec3.jar -O beagle/beagle.21Sep19.ec3.jar'))
    system('chmod +x beagle/beagle.21Sep19.ec3.jar')
    beagle <- "beagle/beagle.21Sep19.ec3.jar"
} else {
    beagle <- "beagle/beagle.21Sep19.ec3.jar"
}



######################################################################
## loadSampleFile
# loadSamplesFile <- function(file, input_folder) {
#     ## debug
#     file = opt$samplesFile; input_folder = opt$inputFolder; #column = opt$samplesColumn
#     ##
#     if (!file.exists(file) ) {
#         write(paste("Sample file",file,"does not exist\n"), stderr())
#         stop()
#     }
#     ### column SAMPLE_ID should be the sample name
#     ### rows can be commented out with #
#     targets <- read.table(file, header = TRUE, as.is = TRUE)
#     for (i in seq.int(nrow(targets$SAMPLE_ID))) {
#         if (targets[i, 1]) {
#             ext <- unique(file_ext(dir(file.path(input_folder,targets[i,1]),pattern = "gz")))
#             if (length(ext) == 0) {
#                 write(paste("Cannot locate fastq or sff file in folder",targets[i,1],"\n"), stderr())
#                 stop()
#             }
#             # targets$type[i] <- paste(ext,sep="/")
#         }
#         else {
#             ext <- file_ext(grep("gz", dir(file.path(input_folder, targets[i, 1])), value = TRUE))
#             if (length(ext) == 0) {
#                 write(paste(targets[i,1],"is not a gz file\n"), stderr())
#                 stop()
#             }
#
#         }
#     }
#     write(paste("samples sheet contains", nrow(targets), "samples to process"),stdout())
#     return(targets)
# }

# #pigz <- system('which pigz 2> /dev/null')
# if (system('which pigz 2> /dev/null', ignore.stdout = TRUE, ignore.stderr = TRUE) == 0) {
#     uncompress <- paste('unpigz', '-p', opt$procs)
# }else{
#     uncompress <- 'gunzip'
# }

######################################################################
## prepareCore
##    Set up the numer of processors to use
##
## Parameters
##    opt_procs: processors given on the option line
##    samples: number of samples
##    targets: number of targets
# if opt_procs set to 0 then expand to samples by targets
prepareCore <- function(opt_procs){
    if (detectCores() < opt$procs) {
        write(paste("number of cores specified (", opt$procs,") is greater than the number of cores available (",detectCores(),")"),stdout())
        paste('Using ', detectCores(), 'threads')
    }
}



######################
# vcfList <- function(samples, reads_folder){
#     mapping_list <- list()
#     for (i in 1:nrow(samples)) {
#         reads <- dir(path = file.path(reads_folder), pattern = "vcf.gz$", full.names = TRUE)
#         map <- lapply(c(".vcf.gz"), grep, x = reads, value = TRUE)
#         names(map) <- c("VCF")
#         map$sampleName <-  samples[i,1]
#         map$VCF <- map$VCF[i]
#         mapping_list[[paste(map$sampleName)]] <- map
#         mapping_list[[paste(map$sampleName, sep = "_")]]
#     }
#     write(paste("Setting up", length(mapping_list), "jobs"),stdout())
#     return(mapping_list)
# }





#samples <- loadSamplesFile(opt$samplesFile, opt$inputFolder)
procs <- prepareCore(opt$procs)
#mapping <- vcfList(samples, opt$inputFolder)



#java -jar ~/programs/beagle.24Aug19.3e8.jar gt=geno_660k_PLINK.converted_from.bed.vcf.gz out=haplotype_0.5Mb_all_lines chrom=5:65603958-66603958



beagle.phase <- function() {
    write(paste('Starting Mapping sample'), stderr())
    try({
        system(paste('java', paste0('-Xmx', opt$memmory,'g'), '-jar', beagle,
                     paste0('gt=', opt$inputFile),
                     paste0('out=', opt$output),
                     paste0('chrom=', opt$window),
                     if (file.exists(opt$map)) {
                         paste0('map=', opt$map)
                     },
                     if (file.exists(opt$excludeSamples)) {
                         opt$excludeSamples
                     },
                     if (file.exists(opt$excludeMarkers)) {
                         opt$excludeMarkers
                     }
        ))})
}
# if (!all(sapply(star.pair.mapping, "==", 0L))) {
#     write(paste("Something went wrong with STAR mapping some jobs failed"),stderr())
#     stop()
# }


run.beagle <- beagle.phase()

# star.pair.mapping <- mclapply(mapping, function(index) {
#     write(paste('Starting Mapping sample', index$sampleName), stderr())
#     try({
#         system(paste('java -jar', beagle,
#                      paste0('gt=', index$VCF),
#                      paste0('out=', index$sampleName, 'output.dat'),
#                      paste0('chrom=',opt$window)
#
#         )})
# }, mc.cores = opt$mprocs
# )
# if (!all(sapply(star.pair.mapping, "==", 0L))) {
#     write(paste("Something went wrong with STAR mapping some jobs failed"),stderr())
#     stop()
# }
