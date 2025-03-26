
library(ggplot2)
library(jsonlite)

exit <- function() { invokeRestart("abort") }

printUsage <- function() {
    cat("lodestarPlots.R v1.0\n");
    cat("----------------------\n\n");
    cat("Written by T. Quinn Smith\n");
    cat("Principal Investigator: Zachary A. Szpiech\n");
    cat("The Pennsylvania State University\n\n");
    cat("Usage: Rscript lodestarPlots.R CMD <windows.json> (<pops.txt>)\n");
    cat("       <windows.json> is the JSON file produced by LODESTAR.\n");
    cat("       <pops.txt> is optional. Labels sample by population.\n");
    cat("CMD:\n");
    cat("----\n");
    cat("mds w i j              Plot component j v. component i of the w'th window.\n");
    cat("axis i                 Plot the values of the i's component along the genome.\n");
    cat("pvals                  Plot the pvalues along the genome.\n");
    cat("tvals                  Plot the t-statistic along the genome.\n");
    cat("\n");
}

axis <- function(windowsJSON, popsFile, i) {
    filename = paste(windowsJSON$output_basename, "_axis_", w, "_", i, "_", j, ".png", sep = "");
    if (popsFile != "") {

    } else {

    }
}

mds <- function(windowsJSON, popsFile, w, i, j) {
    filename = paste(windowsJSON$output_basename, "_mds_", w, "_", i, "_", j, ".png", sep = "");
    if (popsFile != "") {
        points <- as.data.frame(windowsJSON$windows$X[windowsJSON$windows["Window Number"] == w])
        labels <- read.delim(args[2], header = FALSE)[,1];
        df <- data.frame(
            x = points[,1],
            y = points[,2],
            Populations = labels
        );
        p <- ggplot(df, aes(x = x, y = y, color = Populations)) +
            geom_point() + 
            labs(x = paste("MDS Coord ", i), y = paste("MDS Coord", j)) +
            theme(plot.background = element_rect(fill = "white") );
        ggsave(filename, plot = p, width = 6, height = 4, dpi = 300);
    } else {
        df <- as.data.frame(windowsJSON$windows$X[windowsJSON$windows["Window Number"] == w])
        p <- ggplot(df, aes(x = df[,1], y = df[,2])) +
            geom_point() + 
            labs(x = paste("MDS Coord ", i), y = paste("MDS Coord", j)) +
            theme(plot.background = element_rect(fill = "white") );
        ggsave(filename, plot = p, width = 6, height = 4, dpi = 300);
    }
}

cmd <- function(cmd, windowsFile, popsFile, args) {
    windowsJSON = fromJSON(windowsFile);
    switch(cmd,
        mds={
            if (length(args) != 3) {
                cat("mds takes 3 arguments! Exiting!\n");
                return;
            }
            if (args[2] > windowsJSON$k || args[3] > windowsJSON$k || args[2] < 0 || args[3] < 0) {
                cat("Component for MDS plot does not exist! Exiting!\n");
                return;
            }
            if (args[1] > max(unlist(windowsJSON$windows["Window Number"]))) {
                cat("Window for MDS plot does not exist! Exiting!\n");
                return;
            }
            mds(windowsJSON, popsFile, args[1], args[2], args[3]);
        },
        pvals={

        },
        tvals={

        },
        axis={
            if (length(args) != 1) {
                cat("axis only takes one component to graph! Exiting!\n");
                return;
            }
            if (args[1] > windowsJSON$k || args[1] <= 0) {
                cat("axis component is not a valid dimension! Exiting!\n");
                return;
            }
            axis(windowsFile, popsFile, args[1]);
        },
        {

        }
    )
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    printUsage();
    exit();
}

if (args[1] != "mds" && args[1] != "pvals" && args[1] != "tvals" && args[1] != "axis") {
    cat("Unrecognized command! Exiting!\n");
    exit();
}

windowsFile = "";
popsFile = "";

if (substring(args[length(args)], nchar(args[length(args)]) - 3, nchar(args[length(args)])) == ".txt") {
    popsFile = args[length(args)];
    if (!file.exists(popsFile)) {
        cat("<pops.txt> does not exist! Exiting!\n");
        exit();
    }
    windowsFile = args[length(args) - 1];
    if (!file.exists(windowsFile)) {
        cat("<windows.json> does not exist! Exiting!\n");
        exit();
    }
    cmd(args[1], windowsFile, popsFile, args[seq(2, length(args) - 2)]);
} else if (substr(args[length(args)], nchar(args[length(args)]) - 4, nchar(args[length(args)])) == ".json") {
    windowsFile = args[length(args)];
    if (!file.exists(windowsFile)) {
        cat("<windows.json> does not exist! Exiting!\n");
        exit();
    }
    cmd(args[1], windowsFile, popsFile, args[seq(2, length(args) - 1)]);
} else {
    cat("<windows.json> not provided! Exiting!\n");
    exit();
}