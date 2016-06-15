library(minfi)
manifest.file <- "~/Work/cegs/450k/data_files/HumanMethylation450_15017482_v.1.2.csv"
manifest <- minfi:::read.manifest(manifest.file)

## Manifest package

manifestList <- manifest$manifestList
IlluminaHumanMethylation450kmanifest <- do.call(IlluminaMethylationManifest,
                                                list(TypeI = manifestList$TypeI,
                                                     TypeII = manifestList$TypeII,
                                                     TypeControl = manifestList$TypeControl,
                                                     TypeSnpI = manifestList$TypeSnpI,
                                                     TypeSnpII = manifestList$TypeSnpII,
                                                     annotation = "IlluminaHumanMethylation450k"))
save(IlluminaHumanMethylation450kmanifest, # compress = "xz",
     file = "IlluminaHumanMethylation450kmanifest.rda")

## Annotation package

anno <- manifest$manifest
anno$IlmnID <- NULL
nam <- names(anno)
names(nam) <- nam
nam[c("AddressA_ID", "AddressB_ID", "AlleleA_ProbeSeq", "AlleleB_ProbeSeq",
            "Infinium_Design_Type", "Next_Base", "Color_Channel")] <-  c("AddressA", "AddressB",
                                                                         "ProbeSeqA", "ProbeSeqB",
                                                                         "Type", "NextBase", "Color")

names(nam) <- NULL
names(anno) <- nam
rownames(anno) <- anno$Name
anno <- anno[getManifestInfo(IlluminaHumanMethylation450kmanifest, type = "locusNames"),]

Locations <- anno[, c("CHR", "MAPINFO")]
names(Locations) <- c("chr", "pos")
Locations$pos <- as.integer(Locations$pos)
Locations$chr <- paste("chr", Locations$chr, sep = "")
Locations$strand <- ifelse(anno$Strand == "F", "+", "-")
table(Locations$chr, exclude = NULL)
rownames(Locations) <- anno$Name
Locations <- as(Locations, "DataFrame")

Manifest <- anno[, c("Name", "AddressA", "AddressB",
                     "ProbeSeqA", "ProbeSeqB", "Type", "NextBase", "Color")]
Manifest <- as(Manifest, "DataFrame")

SNPs.Illumina <- anno[, c("Probe_SNPs", "Probe_SNPs_10")]
SNPs.Illumina <- as(SNPs.Illumina, "DataFrame")

Islands.UCSC <- anno[, c("UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island")]
names(Islands.UCSC) <- c("Islands_Name", "Relation_to_Island")
Islands.UCSC <- as(Islands.UCSC, "DataFrame")
Islands.UCSC$Relation_to_Island[Islands.UCSC$Relation_to_Island == ""] <- "OpenSea"
table(Islands.UCSC$Relation_to_Island, exclude = NULL)

usedColumns <- c(names(Manifest), names(SNPs.Illumina), 
                 c("CHR", "MAPINFO", "Strand",
                   "Chromosome_36", "Coordinate_36", "Genome_Build"),
                 c("UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island"))
Other <- anno[, setdiff(names(anno), usedColumns)]
Other <- as(Other, "DataFrame")

## We now use an exisitng grSnp object containing a GRanges of relevant SNPs.
## This is created in a separate script

##
## SNP overlap
##

map <- cbind(Locations, Manifest)
map <- GRanges(seqnames = map$chr, ranges = IRanges(start = map$pos, width = 1),
               Strand = map$strand, Type = map$Type)
map <- minfi:::getProbePositionsDetailed(map)
names(map) <- rownames(Locations)

## 137
load("../../../snps/objects/grSnp137CommonSingle.rda")
SNPs.137CommonSingle <- minfi:::.doSnpOverlap(map, grSnp137CommonSingle)
load("../../../snps/objects/grSnp135CommonSingle.rda")
SNPs.135CommonSingle <- minfi:::.doSnpOverlap(map, grSnp135CommonSingle)
load("../../../snps/objects/grSnp132CommonSingle.rda")
SNPs.132CommonSingle <- minfi:::.doSnpOverlap(map, grSnp132CommonSingle)

annoStr <- c(array = "IlluminaHumanMethylation450k",
             annotation = "ilmn12",
             genomeBuild = "hg19")
defaults <- c("Locations", "Manifest", "SNPs.137CommonSingle", "Islands.UCSC", "Other")
annoObj <-
    IlluminaMethylationAnnotation(list(Locations = Locations,
                                       Manifest = Manifest,
                                       SNPs.Illumina = SNPs.Illumina,
                                       SNPs.137CommonSingle = SNPs.137CommonSingle,
                                       SNPs.135CommonSingle = SNPs.135CommonSingle,
                                       SNPs.132CommonSingle = SNPs.132CommonSingle,
                                       Islands.UCSC = Islands.UCSC,
                                       Other = Other
                                       ),
                                  annotation = annoStr, defaults = defaults)
validObject(annoObj)

annoName <- sprintf("%sanno.%s.%s", annoStr["array"], annoStr["annotation"],
                    annoStr["genomeBuild"])
cat("creating object:", annoName, "\n")
assign(annoName, annoObj)
save(list = annoName,
     file = paste(annoName, "rda", sep = "."), compress = "xz")

