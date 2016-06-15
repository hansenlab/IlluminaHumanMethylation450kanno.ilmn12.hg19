processUCSCsnp <- function(snpfile) {
    require(GenomicRanges)
    cat("Reading file\n")
    df <- read.delim(gzfile(snpfile), header = FALSE,
                     stringsAsFactors = FALSE)
    names(df) <- c("chr", "start", "end", "name", "strand",
                   "refNCBI", "class", "alleleFreqs")
    print(table(df$chr))
    cat("Only keeping chrs 1-22, X, Y\n")
    df <- df[df$chr %in% paste0("chr", c(1:22, "X", "Y")),]
    print(table(df$class))
    cat("Only keeping class 'single'\n")
    df <- df[df$class == "single",]
    cat("Computing MAF\n")
    df$alleleFreqs <- sub(",$", "", df$alleleFreqs)
    sp <- strsplit(df$alleleFreqs, ",")
    minFreq <- sapply(sp, function(xx) min(as.numeric(xx)))
    cat("Instatiating object\n")
    grSNP <- GRanges(seqnames = df$chr, strand = df$strand, ranges = IRanges(start = df$start + 1, end = df$end),
                     MAF = minFreq, ref = df$refNCBI)
    names(grSNP) <- df$name
    grSNP
}


grSnp137CommonSingle <- processUCSCsnp("files/snp137Common_small.txt.gz")
save(grSnp137CommonSingle, file = "objects/grSnp137CommonSingle.rda")

grSnp135CommonSingle <- processUCSCsnp("files/snp135Common_small.txt.gz")
save(grSnp135CommonSingle, file = "objects/grSnp135CommonSingle.rda")

grSnp132CommonSingle <- processUCSCsnp("files/snp132Common_small.txt.gz")
save(grSnp132CommonSingle, file = "objects/grSnp132CommonSingle.rda")
