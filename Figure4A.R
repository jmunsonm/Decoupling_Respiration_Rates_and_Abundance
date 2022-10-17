rm(list = ls())

library(plot3D)
library(dplyr)
library(vioplot)

path0 <- "~/Desktop/Git_upload/"
pathDF <- "~/Desktop/Git_upload/Data Files/"
pathFIG <- "~/Desktop/Git_upload/Figures/"
source(paste0(path0,"Figure functions.R"))

#open output plot file
plot_out <- paste0(pathFIG,"Figure_4A.pdf")
pdf(file=plot_out, width=183/25.4, height=10.25)


#set global plotting options/parameters
op <- par(no.readonly = TRUE)
par(mfrow=c(4,2),
    mar=c(3.5,4,1.0,0.5),
    mgp=c(1.75,0.5,0),
    cex.axis=1.5,
    cex.lab=1.75,
    ps=7)

#input color code
input_file <- paste0(pathDF,"Genus_hex_color.csv")
gene_colors <- read.csv(file=input_file, header=TRUE)
gene_colors2 <- gene_colors[1:107,]
gene_colors2 <- rbind(gene_colors2,c("#000000","Unclassified","Unclassified"))
input_file <- paste0(pathDF,"Color_master_key.csv")
class_colors <- read.csv(file=input_file, header=TRUE)
gene_class_colors <- merge(x=gene_colors2, y=class_colors, 
                           by.x="class", by.y="Class", sort=FALSE, all=TRUE)

#rhodopsin vs O2 plot
input_file <- paste0(pathDF,"Fig3F_rhodopsin_per_µm.csv")
Fig3F <- read.csv(file=input_file, header=TRUE)

#assign "class" colors to each point
abc <- merge(x=Fig3F, y=gene_class_colors, by.x="genus", by.y="GTDB_genus", sort=FALSE)

#scale the respiration rate by biovolume
abc$O2_per_volume <- abc$Weighted_avg_O2_consumed_per_cell/abc$spherical.cell.volume..µm3.
abc$count_per_volume <- abc$proteorhodopsin_transcripts_per_cell/abc$spherical.cell.volume..µm3.

#assign plotting points by sample data
abc$shape <- rep(15, nrow(abc))
abc$shape[abc$date == "20190402"] <- 17
abc$shape[abc$date == "20190709"] <- 18

tab_label <- "A"
pt_label <- TRUE


#Figure_4A plot 
#keep only what's needed for plotting
plot_data <- select(abc, "Weighted_avg_O2_consumed_per_cell", "proteorhodopsin_transcripts_per_cell",
                    "ClassColor", "shape", "genus")
colnames(plot_data) <- c("x", "y", "color", "shape", "genus")
#plotting parameters
x_lim <- c(1e-3,1e1)
x_main_label <- expression(paste("Respiration Rate (fmol O"[2]," cell"^-1," h"^-1*")"))
major_x_tic_location <- c(1e-3,1e-2,1e-1,1e0,1e1)
major_x_tic_labels <- c("0.001", "0.01", "0.1", "1", "10")
minor_x_tic_location <- c((2:9)*1e-3, (2:9)*1e-2, (2:9)*1e-1, (2:9)*1e0)

y_lim <- c(1e-3,1e1)
y_main_label <- expression(paste("Rhodopsin transcripts (cell"^-1*")"))
major_y_tic_location <- 10^(-3:1)
major_y_tic_labels <- c("0.001 ","0.01 ", "0.1 ","1 ", "10 ")
minor_y_tic_location <- c((2:9)*1e-3, (2:9)*1e-2, (2:9)*1e-1, (2:9)*1e0)

reg_eqn_xy <- NULL
label_points_xy <- rbind(c(0.02, 0.02, 0.01, 0.4),c(1.00,0.03, 3.00, 0.02))
tab_y <- 10^(1+0.075*(4))
par(mfg=c(1,1))

#label_points_xy <- rbind(c(0.01, 0.03, 0.1, 0.4),c(9, 7.5, 5.5, 4.5))

plot_logxy(plot_data,
           x_lim, x_main_label, major_x_tic_location, major_x_tic_labels, minor_x_tic_location,
           y_lim, y_main_label, major_y_tic_location, major_y_tic_labels, minor_y_tic_location,
           reg_eqn_xy, pt_label, label_points_xy, tab_y, tab_label)

par(op)
dev.off()

