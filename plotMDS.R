
# File: plotPeaks.R 
# Date: 14 January 2025
# Author: T. Quinn Smith
# Principal Investigator: Dr. Zachary A. Szpiech
# Purpose: Plot points in 2 dimensions from a tsv file of points.

# Our only dependency.
library(ggplot2)

main <- function() {
    # Get command line arguments.
    args <- commandArgs(trailingOnly = TRUE)
    # If a filename and/or -p option was not given, print error and exit.
    if (length(args) != 1 && length(args) != 2) {
        cat("Usage: Rscript plotMDS.R [-p] <summary.tsv>\n")
    } else {
        fileName = args[1]
        if (!file.exists(fileName)) {
            cat("The file", fileName, "does not exist!\n")
        } else {
            # Create our plot.
            df <- read.delim(fileName, header = FALSE)
            p <- ggplot(df, aes(x = df[,1], y = df[,2])) +
                geom_point() + 
                labs(x = "MDS Coord 1", y = "MDS Coord 2") +
                theme(plot.background = element_rect(fill = "white") )
            ggsave(paste(substr(fileName, 1, nchar(fileName) - 4), ".png", sep = ""), plot = p, width = 6, height = 4, dpi = 300)
        }
    }
}

main()