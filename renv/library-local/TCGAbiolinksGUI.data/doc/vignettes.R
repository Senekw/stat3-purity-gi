## ----eval = FALSE-------------------------------------------------------------
# if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")
# BiocManager::install(("TCGAbiolinksGUI.data")

## ----eval = FALSE-------------------------------------------------------------
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")

## ----eval = TRUE--------------------------------------------------------------
library(TCGAbiolinksGUI.data)

## ----eval=FALSE, include=TRUE-------------------------------------------------
# # Defining parameters
# getGDCdisease <- function(){
#   projects <- TCGAbiolinks:::getGDCprojects()
#   projects <- projects[projects$id != "FM-AD",]
#   disease <-  projects$project_id
#   idx <- grep("disease_type",colnames(projects))
#   names(disease) <-  paste0(projects[[idx]], " (",disease,")")
#   disease <- disease[sort(names(disease))]
#   return(disease)
# }

## -----------------------------------------------------------------------------
data(GDCdisease)
DT::datatable(as.data.frame(GDCdisease))

## ----eval=FALSE, include=TRUE-------------------------------------------------
# getMafTumors <- function(){
#   root <- "https://gdc-api.nci.nih.gov/data/"
#   maf <- fread("https://gdc-docs.nci.nih.gov/Data/Release_Notes/Manifests/GDC_open_MAFs_manifest.txt",
#                data.table = FALSE, verbose = FALSE, showProgress = FALSE)
#   tumor <- unlist(lapply(maf$filename, function(x){unlist(str_split(x,"\\."))[2]}))
#   proj <- TCGAbiolinks:::getGDCprojects()
# 
#   disease <-  gsub("TCGA-","",proj$project_id)
#   idx <- grep("disease_type",colnames(proj))
#   names(disease) <-  paste0(proj[[idx]], " (",proj$project_id,")")
#   disease <- sort(disease)
#   ret <- disease[disease %in% tumor]
#   return(ret)
# }

## -----------------------------------------------------------------------------
data(maf.tumor)
DT::datatable(as.data.frame(maf.tumor))

## ----eval = FALSE, include = TRUE---------------------------------------------
# library(readr)
# library(readxl)
# library(dplyr)
# library(caret)
# library(randomForest)
# library(doMC)
# library(e1071)
# 
# # Control random number generation
# set.seed(210) # set a seed to RNG
# 
# # register number of cores to be used for parallel evaluation
# registerDoMC(cores = parallel::detectCores())

## ----eval=FALSE, include=TRUE-------------------------------------------------
# file <- "https://tcga-data.nci.nih.gov/docs/publications/lgggbm_2016/LGG.GBM.meth.txt"
# if(!file.exists(basename(file))) downloader::download(file,basename(file))
# LGG.GBM <- as.data.frame(readr::read_tsv(basename(file)))
# rownames(LGG.GBM) <- LGG.GBM$Composite.Element.REF
# idx <- grep("TCGA",colnames(LGG.GBM))
# colnames(LGG.GBM)[idx] <- substr(colnames(LGG.GBM)[idx], 1, 12) # reduce complete barcode to sample identifier (first 12 characters)

## ----eval=FALSE, include=TRUE-------------------------------------------------
# file <- "http://www.cell.com/cms/attachment/2045372863/2056783242/mmc2.xlsx"
# if(!file.exists(basename(file))) downloader::download(file,basename(file))
# metadata <-  readxl::read_excel(basename(file), sheet = "S1A. TCGA discovery dataset", skip = 1)
# DT::datatable(
#   metadata[,c("Case",
#               "Pan-Glioma DNA Methylation Cluster",
#               "Supervised DNA Methylation Cluster",
#               "IDH-specific DNA Methylation Cluster")]
# )

## ----eval=FALSE, include=TRUE-------------------------------------------------
# file <- "http://zwdzwd.io/InfiniumAnnotation/current/EPIC/EPIC.manifest.hg38.rda"
# if(!file.exists(basename(file))) downloader::download(file,basename(file))
# load(basename(file))

## ----eval=FALSE, include=TRUE-------------------------------------------------
# file <- "https://tcga-data.nci.nih.gov/docs/publications/lgggbm_2015/PanGlioma_MethylationSignatures.xlsx"
# if(!file.exists(basename(file))) downloader::download(file,basename(file))

## ----eval=FALSE, include=TRUE-------------------------------------------------
# sheet <- "1,300 pan-glioma tumor specific"
# trainingset <- grep("mut|wt",unique(metadata$`Pan-Glioma DNA Methylation Cluster`),value = T)
# trainingcol <- c("Pan-Glioma DNA Methylation Cluster")

## ----eval=FALSE, include=TRUE-------------------------------------------------
# plat <- "EPIC"
# signature.probes <-  read_excel("PanGlioma_MethylationSignatures.xlsx",  sheet = sheet)  %>% pull(1)
# samples <- dplyr::filter(metadata, `IDH-specific DNA Methylation Cluster` %in% trainingset)
# RFtrain <- LGG.GBM[signature.probes, colnames(LGG.GBM) %in% as.character(samples$Case)] %>% na.omit

## ----eval=FALSE, include=TRUE-------------------------------------------------
# RFtrain <- RFtrain[!EPIC.manifest.hg38[rownames(RFtrain)]$MASK.general,]

## ----eval=FALSE, include=TRUE-------------------------------------------------
# trainingdata <- t(RFtrain)
# trainingdata <- merge(trainingdata, metadata[,c("Case", trainingcol[model])], by.x=0,by.y="Case", all.x=T)
# rownames(trainingdata) <- as.character(trainingdata$Row.names)
# trainingdata$Row.names <- NULL

## ----eval=FALSE, include=TRUE-------------------------------------------------
# nfeat <- ncol(trainingdata)
# trainingdata[,trainingcol] <-  factor(trainingdata[,trainingcol])
# mtryVals <- floor(sqrt(nfeat))
# for(i in floor(seq(sqrt(nfeat), nfeat/2, by = 2 * sqrt(nfeat)))) {
#   print(i)
#   x <- as.data.frame(
#     tuneRF(
#       trainingdata[,-grep(trainingcol[model],colnames(trainingdata))],
#       trainingdata[,trainingcol[model]],
#       stepFactor=2,
#       plot= FALSE,
#       mtryStart = i
#     )
#   )
#   mtryVals <- unique(c(mtryVals, x$mtry[which (x$OOBError == min(x$OOBError))]))
# }
# mtryGrid <- data.frame(.mtry = mtryVals)

## ----eval=FALSE, include=TRUE-------------------------------------------------
# fitControl <- trainControl(
#   method = "repeatedcv",
#   number = 10,
#   verboseIter = TRUE,
#   repeats = 10
# )

## ----eval=FALSE, include=TRUE-------------------------------------------------
# glioma.idh.model <- train(
#   y = trainingdata[,trainingcol], # variable to be trained on
#   x = trainingdata[,-grep(trainingcol,colnames(trainingdata))], # Daat labels
#   data = trainingdata, # Data we are using
#   method = "rf", # Method we are using
#   trControl = fitControl, # How we validate
#   ntree = 5000, # number of trees
#   importance = TRUE,
#   tuneGrid = mtryGrid, # set mtrys, the value that procuded a better model is used
# )

## -----------------------------------------------------------------------------
data(glioma.idh.model)
glioma.idh.model

## ----eval=FALSE, include=TRUE-------------------------------------------------
# sheet <- "1,308 IDHmutant tumor specific "
# trainingset <- grep("mut",unique(metadata$`IDH-specific DNA Methylation Cluster`),value = T)
# trainingcol <- "IDH-specific DNA Methylation Cluster"

## -----------------------------------------------------------------------------
data(glioma.idhmut.model)
glioma.idhmut.model

## ----eval=FALSE, include=TRUE-------------------------------------------------
# sheet <- "163  probes that define each TC"
# trainingset <- c("IDHmut-K1","IDHmut-K2")
# trainingcol <- "Supervised DNA Methylation Cluster"

## -----------------------------------------------------------------------------
data("glioma.gcimp.model")
glioma.gcimp.model

## ----eval=FALSE, include=TRUE-------------------------------------------------
# sheet <- "914 IDHwildtype tumor specific "
# trainingset <- grep("wt",unique(metadata$`IDH-specific DNA Methylation Cluster`),value = T))
# trainingcol <- "IDH-specific DNA Methylation Cluster"

## -----------------------------------------------------------------------------
data("glioma.idhwt.model")
glioma.idhwt.model

## -----------------------------------------------------------------------------
data("probes2rm")
head(probes2rm)

## ----eval = FALSE-------------------------------------------------------------
# scraplinks <- function(url){
#     # Create an html document from the url
#     webpage <- xml2::read_html(url)
#     # Extract the URLs
#     url_ <- webpage %>%
#         rvest::html_nodes("a") %>%
#         rvest::html_attr("href")
#     # Extract the link text
#     link_ <- webpage %>%
#         rvest::html_nodes("a") %>%
#         rvest::html_text()
#     tb <- tibble::tibble(link = link_, url = url_)
#     tb <- tb %>% dplyr::filter(tb$link == "Download")
#     return(tb)
# }
# 
# library(htmltab)
# library(dplyr)
# library(tidyr)
# root <- "http://linkedomics.org"
# root.download <- file.path(root,"data_download")
# linkedOmics <- htmltab(paste0(root,"/login.php#dataSource"))
# linkedOmics <- linkedOmics %>% unite(col = "download_page","Cohort Source","Cancer ID", remove = FALSE,sep = "-")
# linkedOmics.data <- plyr::adply(linkedOmics$download_page,1,function(project){
#     url <- file.path(root.download,project)
#     tryCatch({
#         tb <- cbind(tibble::tibble(project),htmltab(url),scraplinks(url))
#         tb$Link <- tb$link <- NULL
#         tb$dataset <- gsub(" \\(.*","",tb$`OMICS Dataset`)
#         tb
#     }, error = function(e) {
#         message(e)
#         return(NULL)
#     }
#     )
# },.progress = "time",.id = NULL)
# usethis::use_data(linkedOmics.data,internal = FALSE,compress = "xz")

## ----eval = FALSE-------------------------------------------------------------
# get_gene_information_biomart <- function(
#     genome = c("hg38","hg19"),
#     TSS = FALSE
# ){
#     requireNamespace("biomaRt")
#     genome <- match.arg(genome)
#     tries <- 0L
#     msg <- character()
#     while (tries < 3L) {
#         gene.location <- tryCatch({
#             host <- ifelse(
#                 genome == "hg19",
#                 "grch37.ensembl.org",
#                 "www.ensembl.org"
#             )
#             mirror <- list(NULL, "useast", "uswest", "asia")[[tries + 1]]
#             ensembl <- tryCatch({
#                 message(
#                     ifelse(
#                         is.null(mirror),
#                         paste0("Accessing ", host, " to get gene information"),
#                         paste0("Accessing ", host, " (mirror ", mirror, ")")
#                     )
#                 )
#                 biomaRt::useEnsembl(
#                     "ensembl",
#                     dataset = "hsapiens_gene_ensembl",
#                     host = host,
#                     mirror = mirror
#                 )
#             }, error = function(e) {
#                 message(e)
#                 return(NULL)
#             })
# 
#             # Column values we will recover from the database
#             attributes <- c(
#                 "ensembl_gene_id",
#                 "external_gene_name",
#                 "entrezgene",
#                 "chromosome_name",
#                 "strand",
#                 "end_position",
#                 "start_position",
#                 "gene_biotype"
#             )
# 
#             if (TSS)  attributes <- c(attributes, "transcription_start_site")
# 
#             db.datasets <- biomaRt::listDatasets(ensembl)
#             description <- db.datasets[db.datasets$dataset == "hsapiens_gene_ensembl", ]$description
#             message(paste0("Downloading genome information (try:", tries, ") Using: ", description))
#             gene.location <- biomaRt::getBM(
#                 attributes = attributes,
#                 filters = "chromosome_name",
#                 values = c(seq_len(22),"X","Y"),
#                 mart = ensembl
#             )
#             gene.location
#         }, error = function(e) {
#             msg <<- conditionMessage(e)
#             tries <<- tries + 1L
#             NULL
#         })
#         if (!is.null(gene.location)) break
#         if (tries == 3L) stop("failed to get URL after 3 tries:", "\n  error: ", msg)
#     }
# }
# gene.location.hg19 <- get_gene_information_biomart("hg19")
# save(gene.location.hg19, file = "gene.location.hg19.rda")
# 
# gene.location.hg38 <- get_gene_information_biomart("hg38")
# save(gene.location.hg38, file = "gene.location.hg38.rda")

## ----eval = FALSE-------------------------------------------------------------
# library(biomaRt)
# getTSS <- function(
#     genome = "hg38",
#     TSS = list(upstream = NULL, downstream = NULL)
# ) {
#     host <- ifelse(genome == "hg19",  "grch37.ensembl.org", "www.ensembl.org")
#     ensembl <- tryCatch({
#         useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host =  host)
#     },  error = function(e) {
#         useEnsembl(
#           biomart = "ensembl",
#             dataset = "hsapiens_gene_ensembl",
#             host =  host
#         )
#     })
#     attributes <- c(
#         "chromosome_name",
#         "start_position",
#         "end_position",
#         "strand",
#         "transcript_start",
#         "transcription_start_site",
#         "transcript_end",
#         "ensembl_transcript_id",
#         "ensembl_gene_id",
#         "external_gene_name"
#     )
# 
#     chrom <- c(1:22, "X", "Y")
#     db.datasets <- listDatasets(ensembl)
#     description <- db.datasets[db.datasets$dataset == "hsapiens_gene_ensembl", ]$description
#     message(paste0("Downloading transcripts information. Using: ", description))
# 
#     tss <- getBM(
#         attributes = attributes,
#         filters = c("chromosome_name"),
#         values = list(chrom),
#         mart = ensembl
#     )
#     tss <- tss[!duplicated(tss$ensembl_transcript_id), ]
#     if (genome == "hg19") tss$external_gene_name <- tss$external_gene_id
#     tss$chromosome_name <-  paste0("chr", tss$chromosome_name)
#     tss$strand[tss$strand == 1] <- "+"
#     tss$strand[tss$strand == -1] <- "-"
# 
#     tss <- makeGRangesFromDataFrame(
#         tss,
#         start.field = "transcript_start",
#         end.field = "transcript_end",
#         keep.extra.columns = TRUE
#     )
# 
#     if (!is.null(TSS$upstream) & !is.null(TSS$downstream)) {
#         tss <- promoters(
#             tss,
#             upstream = TSS$upstream,
#             downstream = TSS$downstream
#         )
#     }
#     return(tss)
# }
# 
# gene.location.hg19.tss <- getTSS("hg19")
# save(gene.location.hg19.tss, file = "gene.location.hg19.tss.rda")
# 
# gene.location.hg38.tss <- getTSS("hg38")
# save(gene.location.hg38.tss, file = "gene.location.hg38.tss.rda")
# 

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

