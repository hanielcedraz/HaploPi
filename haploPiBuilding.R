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
#suppressPackageStartupMessages(library("dplyr"))

option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = "file.txt",
                help = "The filename of the sample file [default %default]",
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
    make_option(c("-b", "--haplotypeFolder"), type = "character", default = '01-haplotypeFolder',
                help = "Directory where to store the mapping results [default %default]",
                dest = "haplotypeFolder"),
    make_option(c('-s', '--snp'), type  =  'character', default  =  "FALSE",
                help = 'SNP-ID to filter haplotypes [default %default]',
                dest = 'snpID'),
    make_option(c("-w", "--window"), type = "integer", default = 30,
                help = "values to left and right of SNP position",
                dest = "window"),
    # make_option(c("-g", "--gtfTargets"), type = "character", default = "gtf_targets.gtf",
    #             help = "Path to a gtf file, or tab delimeted file with [target gtf] to run mapping against. If would like to run without gtf file, -g option is not required [default %default]",
    #             dest = "gtfTarget"),
    make_option(c("-p", "--processors"), type = "integer", default = 8,
                help = "number of processors to use [defaults %default]",
                dest = "procs"),
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


# if (!file.exists('beagle')) {
#     dir.create(file.path('beagle'), recursive = TRUE, showWarnings = FALSE)
#     system(paste('wget https://faculty.washington.edu/browning/beagle/beagle.21Sep19.ec3.jar -O beagle/beagle.21Sep19.ec3.jar'))
# }





prepareCore <- function(opt_procs){
    if (detectCores() < opt$procs) {
        write(paste("number of cores specified (", opt$procs,") is greater than the number of cores available (",detectCores(),")"),stdout())
        paste('Using ', detectCores(), 'threads')
    }
}

procs <- prepareCore(opt$procs)




input <- opt$inputFile
phenotype <- opt$phenotypeFile
output1 <- opt$output
outputprocessed <- paste0('processed_', opt$output)
final_output <- paste0('built_', opt$output)


# input <- 'haplotype_phased_LLLL_chr5.vcf.gz'
# phenotype <- 'corrected_phenotypes_LLLL_backfat.txt'
# output1 <- 'haplotype_phased_built_LLLL_chr5.dat'
# outputprocessed <- paste0('processed_', input)
# final_output <- paste0('built_', output1)
cat('\n')
write(paste('Loading datasets'), stderr())
cat('\n')

haplotype <- fread(input, header = TRUE, check.names = FALSE, sep = "\t", drop = c(1, 2, 3, 6:9), nThread = opt$procs)
haplotype[1:10,1:10]


snp <- fread(input, header = TRUE, check.names = FALSE, sep = "\t", select = c(3, 2), nThread = opt$procs)


cat('\n')
write(paste('saving output'), stderr())
cat('\n')
fwrite(haplotype, output1, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t", nThread = opt$procs)

cat('\n')
write(paste('processing text. replacing underscore and pipe to tab using tr function'), stderr())
cat('\n')
system(paste("cat", output1, "| tr '_' '\t' | tr '|' '\t' >", outputprocessed))


cat('\n')
write(paste('converting haplotypes 0 and 1 to REF and ALT'), stderr())
cat('\n')

haplo <- as.matrix(fread(outputprocessed, header = TRUE, check.names = FALSE, sep = "\t", fill = TRUE), nThread = opt$procs)

cat('\n')
write(paste('haplotype before convert'), stderr())
haplo[1:10,1:10]
cat('\n')





for (i in 1:nrow(haplo)) {
    for (j in 3:ncol(haplo)) {
        haplo[i,j][haplo[i,j] == 0] <- haplo[i,1]
        haplo[i,j][haplo[i,j] == 1] <- haplo[i,2]
        # haplo[i,j][haplo[i,j] == 0] <- haplo[i,1]
        # haplo[i,j][haplo[i,j] == 1] <- haplo[i,2]
    }
}

cat('\n')
write(paste('haplotype after convert'), stderr())
haplo[1:10,1:10]
cat('\n')



cat('\n')
write(paste('Creating final datasets'), stderr())
cat('\n')
final_df <- data.table(SNP = snp, haplo, keep.rownames = TRUE)
final_df[1:10,1:10]



if (opt$snpID != 'FALSE') {
    leadsnp <- which(grepl(paste(opt$snpID), final_df$SNP.ID))
    minor <- leadsnp - opt$window
    major <- leadsnp + opt$window
    final_df <- final_df[minor:major,]
}
# leadsnp <- which(grepl("AX-116265245", final_df$SNP.ID))
# minor <- leadsnp - 20
# major <- leadsnp + 20
# final_df <- final_df[minor:major,]
#AX-116265245


cat('\n')
transfinaldf <- transpose(final_df)
transfinaldf[1:10,1:10]
cat('\n')

pheno <- fread(phenotype, header = FALSE, check.names = FALSE, sep = "\t", select = c(1, 4), dec = ',', nThread = opt$procs, col.names = c('ID', 'phenotype'))
pheno <- pheno[rep(1:nrow(pheno),each = 2),]
extradata <- data.table(ID = c('SNP.ID', 'SNP.POS', 'REF', 'ALT'), phenotype = c('NA', 'NA', 'NA', 'NA'))

pheno <- rbind(extradata, pheno)
cat('\n')

finalDF <- data.table(pheno[,1], transfinaldf, pheno[,2], keep.rownames = TRUE)

cat('\n')
finalDF[1:10,c(1:10)]
cat('\n')


cat('\n')
write(paste('Saving final haplotype'), stderr())
cat('\n')
fwrite(finalDF, paste0('transposed_', final_output,'.csv'), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t", nThread = opt$procs)


if (file.exists(paste0('transposed_', final_output,'.csv'))) {
    unlink(c(output1, outputprocessed), recursive = TRUE)
}
#write.zoo(finalDF, paste0('transposed_', final_output,'.csv'), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

ncoll <- ncol(finalDF) - 1
nrown <- nrow(finalDF) - 4


uni <- data.table(finalDF[5:nrow(finalDF),2:ncoll])



unihaplo <- unite(uni, col = 'Haplotypes', sep = '')
unihaplo



df <- data.table(finalDF[5:nrow(finalDF),1], unihaplo, pheno[5:nrow(pheno),2])
counting <- count(df, vars = Haplotypes, sort = TRUE, name = 'Number_of_haplotypes')
freq <- list()
for (i in 1:nrow(counting)) {
    freq[i] <- (counting[i,2]/nrown)*100
}
freq <- data.table(freq)
counting <- data.table(counting, freq)
colnames(counting) <- c('Haplotypes', 'Number', 'freq')

#frequency = totalCount/totalpopulation

fwrite(df, paste0('join_haplotype_', final_output, '.csv'), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t", nThread = opt$procs)

fwrite(counting, paste0('haplotype_frequency_', final_output, '.csv'), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t", nThread = opt$procs)


