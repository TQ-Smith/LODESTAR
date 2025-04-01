
# File: lodestarPlots.R
# Date: 27 March 2025
# Author: T. Quinn Smith
# Principal Investigator: Dr. Zachary A. Szpiech
# Purpose: Create plots from LODESTAR analysis.

library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
library(jsonlite)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)

exit <- function() { invokeRestart("abort") }

# Print the usage statement for the script.
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
    cat("tvals                  Plot the t-statistic along the genome. Ignores <pops.txt>.\n");
    cat("print w                Prints the coordinates of the w'th window.\n")
    cat("\n");
}

# Plot the t-statistic along the genome.
tvals <- function(windowsJSON) {
    # Drop invalid windows.
    windowsJSON$windows = windowsJSON$windows[windowsJSON$windows["Chromosome"] != "Global" & windowsJSON$windows["t-statistic"] != -1,];
    filename = paste(windowsJSON$output_basename, "_tvals.png", sep = "");
    ylab = ""
    if (windowsJSON$similarity) {
        ylab = expression("t"["s"]);
    } else { 
        ylab = expression("t"["d"]);
    }
    data <- data.frame(
        CHR = windowsJSON$windows["Chromosome"],  
        BP = as.numeric(windowsJSON$windows["Start Coordinate"][,1]),  
        T = as.numeric(windowsJSON$windows["t-statistic"][,1]),
        SNP = "."
    );
    colnames(data) <- c("CHR", "BP", "T", "SNP");
    data <- data %>%
        mutate(CHR = as.numeric(gsub("chr", "", CHR)));
    # https://r-graph-gallery.com/101_Manhattan_plot.html
    don <- data %>% 
        group_by(CHR) %>% 
        summarise(chr_len=max(BP)) %>% 
        mutate(tot=cumsum(chr_len)-chr_len) %>%
        select(-chr_len) %>%
        left_join(data, ., by=c("CHR"="CHR")) %>%
        arrange(CHR, BP) %>%
        mutate( BPcum=BP+tot);
    axisdf = don %>%
        group_by(CHR) %>%
        summarize(center=( max(BPcum) + min(BPcum) ) / 2 );

    plot <- ggplot(don, aes(x=BPcum, y=T)) +
        geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
        scale_color_manual(values = rep(c("red", "blue"), 22 )) +
        scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
        scale_y_continuous(expand = c(0, 0)) + 
        theme_bw() +
        theme( 
            legend.position="none",
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()
        ) +
        labs(
            title = windowsJSON$input_file,
            x = "Chromosome Position",
            y = ylab,
        ) + ylim(0, 1);
    ggsave(filename, plot, width = 10, height = 6, dpi = 300);
}

# Plot an axis's value at each window for all samples.
axis <- function(windowsJSON, popsFile, i) {
    windowsJSON$windows = windowsJSON$windows[windowsJSON$windows["t-statistic"] != -1,];
    filename = paste(windowsJSON$output_basename, "_axis_", i, ".png", sep = "");

    # Expand chrom and position for all samples in every window.
    numSamples = length(windowsJSON$samples);
    chroms = unlist(lapply(windowsJSON$windows$Chromosome, function(x) rep(x, numSamples)));
    bp = unlist(lapply(windowsJSON$windows["Start Coordinate"], function(x) rep(x, numSamples)));
    components = unlist(lapply(windowsJSON$windows$X, function(x) x[,as.integer(i)]));

    # If we have population membership, then label the points.
    if (popsFile != "") {
        labels <- read.delim(popsFile, header = FALSE)[,1];
        labels <- unlist(rep(labels, length(windowsJSON$windows$X)));
        data <- data.frame(
            CHR = chroms,
            BP = bp,
            COMP = components,
            Populations = labels
        );
        g <- ggplot(data, aes(x = BP, y = COMP, color = as.factor(Populations))) +
            geom_point( aes(color=as.factor(Populations)), alpha=0.8, size=1.3) +
            labs(x = "Chromosome", y = paste("Axis ", i), title=windowsJSON$input_file) +
            theme_bw() +
            theme( 
                legend.position="none",
                panel.border = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank()
            )
        ggsave(filename=filename, plot=g, width = 10, height = 5, dpi = 300);
    } else {
        data <- data.frame(
            CHR = chroms,
            BP = bp,
            COMP = components
        );
        g <- ggplot(data, aes(x = BP, y = COMP)) +
            geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
            scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
            labs(x = "Chromosome", y = paste("Axis ", i), title=windowsJSON$input_file) +
            theme_bw() +
            theme( 
                legend.position="none",
                panel.border = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank()
            )
        ggsave(filename=filename, plot=g, width = 10, height = 5, dpi = 300);
    }
}

# Plot MDS for the i and j axes for the given window.
mds <- function(windowsJSON, popsFile, w, i, j) {
    filename = paste(windowsJSON$output_basename, "_mds_", w, "_", i, "_", j, ".png", sep = "");
    if (popsFile != "") {
        points <- as.data.frame(windowsJSON$windows$X[windowsJSON$windows["Window Number"] == w])
        labels <- read.delim(popsFile, header = FALSE)[,1];
        df <- data.frame(x = points[,as.integer(i)], y = points[,as.integer(j)], Populations = labels);
        p <- ggplot(df, aes(x = x, y = y, color = Populations)) +
            geom_point() + 
            labs(x = paste("Axis ", i), y = paste("Axis ", j)) +
            theme(plot.background = element_rect(fill = "white") );
        ggsave(filename, plot = p, width = 6, height = 4, dpi = 300);
    } else {
        points <- as.data.frame(windowsJSON$windows$X[windowsJSON$windows["Window Number"] == w])
        df <- data.frame(x = points[,as.integer(i)], y = points[,as.integer(j)]);
        p <- ggplot(df, aes(x = x, y = y)) +
            geom_point() + 
            labs(x = paste("Axis ", i), y = paste("Axis ", j)) +
            theme(plot.background = element_rect(fill = "white") );
        ggsave(filename, plot = p, width = 6, height = 4, dpi = 300);
    }
}

# Print the coordinates of the window.
printWindow <- function(windowsJSON, w) {
    points <- as.data.frame(windowsJSON$windows$X[windowsJSON$windows["Window Number"] == w]);
    write.table(points, sep = "\t", row.names = FALSE, col.names = FALSE);
}

# Warning messages and run a give command.
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
            if (args[1] != 0 && windowsJSON$windows["t-statistic"][windowsJSON$windows["Window Number"] == args[1]] == -1) {
                cat("Window for MDS plot does not exist! Exiting!\n");
                return;
            }
            mds(windowsJSON, popsFile, args[1], args[2], args[3]);
        },
        tvals={
            tvals(windowsJSON);
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
            axis(windowsJSON, popsFile, args[1]);
        },
        print={
            if (args[1] != 0 && windowsJSON$windows["t-statistic"][windowsJSON$windows["Window Number"] == args[1]] == -1) {
                cat("Window for MDS plot does not exist! Exiting!\n");
                return;
            }
            printWindow(windowsJSON, args[1]);
        },
        {
            return;
        }
    )
}

# Parse input arguments.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    printUsage();
    exit();
}

if (args[1] != "print" && args[1] != "mds" && args[1] != "tvals" && args[1] != "axis") {
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