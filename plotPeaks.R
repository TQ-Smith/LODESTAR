
# File: plotPeaks.R 
# Date: 14 January 2025
# Author: T. Quinn Smith
# Principal Investigator: Dr. Zachary A. Szpiech
# Purpose: Plot the results of a LODESTAR analysis in a Manhattan plot like fashion.

# Our only dependency.
library(qqman)

main <- function() {
    # Get command line arguments.
    args <- commandArgs(trailingOnly = TRUE)
    # If a filename and/or -p option was not given, print error and exit.
    if (length(args) != 1 && length(args) != 2) {
        cat("Usage: Rscript plotPeaks.R [-p] <summary.tsv>\n")
    } else {
        # Get what we should plot (t or p) and file name from user.
        useT = TRUE
        fileName = args[1]
        if (args[1] == "-p") {
            useT = FALSE
            fileName = args[2]
        }
        if (!file.exists(fileName)) {
            cat("The file", fileName, "does not exist!\n")
        } else {
            # Create our plot.
            df <- read.delim(fileName, skip=1, header=TRUE)
            df <- df[!df$t.stat == -1,]
            df[9] <- lapply(df[9], function (x) 1 - x)
            manhat <- data.frame(
                SNP <- ".",
                CHR <- as.numeric(factor(df$Chr, levels = unique(df$Chr))),
                BP <- df["Start"],
                P <- df["t.stat"]
            )
            if (useT) {
                # Take difference to observe magnitude.
                print(paste("Output file: ", substr(fileName, 1, nchar(fileName) - 4), "_t.png", sep = ""))
            } else {
                manhat$P <- df["p.val"]
                print(paste("Output file: ", substr(fileName, 1, nchar(fileName) - 4), "_p.png", sep = ""))
            }
            colnames(manhat) <- c("SNP", "CHR", "BP", "P")
            # Save plot to png.
            if (useT && length(unique(manhat$CHR)) == 1) {
                # Single chromosome.
                png(paste(substr(fileName, 1, nchar(fileName) - 4), "_t.png", sep = ""), width = 1600, height = 1200, res = 300)
                m <- manhattan(manhat, chr = "CHR", bp = "BP", p = "P", snp = "SNP", ylab = "-log(1-t)", main = substr(fileName, 1, length(fileName) - 4), xlim = c(min(manhat$BP) - 5, max(manhat$BP) + 1), ylim = c(0, max(manhat$P) + 2))
            } else if (!useT && length(unique(manhat$CHR)) == 1) {
                # Single chromosome.
                png(paste(substr(fileName, 1, nchar(fileName) - 4), "_p.png", sep = ""), width = 1600, height = 1200, res = 300)
                m <- manhattan(manhat, chr = "CHR", bp = "BP", p = "P", snp = "SNP", main = substr(fileName, 1, length(fileName) - 4), xlim = c(min(manhat$BP) - 5, max(manhat$BP) + 1), ylim = c(0, max(manhat$P) + 2))
            } else if (useT) {
                png(paste(substr(fileName, 1, nchar(fileName) - 4), "_t.png", sep = ""), width = 1600, height = 1200, res = 300)
                m <- manhattan(manhat, chr = "CHR", bp = "BP", p = "P", snp = "SNP", ylab = "-log(1-t)", main = substr(fileName, 1, length(fileName) - 4))
            } else {
                png(paste(substr(fileName, 1, nchar(fileName) - 4), "_p.png", sep = ""), width = 1600, height = 1200, res = 300)
                m <- manhattan(manhat, chr = "CHR", bp = "BP", p = "P", snp = "SNP", main = substr(fileName, 1, length(fileName) - 4))
            }
            dev.off()
        }
    }
}

main()