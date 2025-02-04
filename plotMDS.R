
# File: plotMDS.R 
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
        cat("Usage: Rscript plotMDS.R <coords.tsv> (<pops.txt>)\n")
    } else {
        fileName = args[1]
        if (!file.exists(fileName)) {
            cat("The file", fileName, "does not exist!\n")
        } else {
            # If populations are specified.
            if (length(args) == 2) {
                if (!file.exists(args[2])) {
                    cat("<pops.txt> does not exist!\n")
                } else {
                    points <- read.delim(fileName, header = FALSE)
                    labels <- read.delim(args[2], header = FALSE)[,1]
                    df <- data.frame(
                        x = points[,1],
                        y = points[,2],
                        Populations = labels
                    )
                    p <- ggplot(df, aes(x = x, y = y, color = Populations)) +
                        geom_point() + 
                        labs(x = "MDS Coord 1", y = "MDS Coord 2") +
                        theme(plot.background = element_rect(fill = "white") )
                    ggsave(paste(substr(fileName, 1, nchar(fileName) - 4), ".png", sep = ""), plot = p, width = 6, height = 4, dpi = 300)
                }
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
}

main()