rm(list = ls())

library(plot3D)
library(dplyr)
library(vioplot)

path0 <- "~/Desktop/Git_upload/"
pathDF <- "~/Desktop/Git_upload/Data Files/"
pathFIG <- "~/Desktop/Git_upload/Figures/"
source(paste0(path0,"Figure functions.R"))

pt_label <- TRUE   #use TRUE if you want selected points labeled, FALSE otherwise
min_O2_keep <- 0.000
min_proteorhodopsin <- 0
max_proteorhodopsin <- 10

x_adj <- 0.50   #location of regression equation (% distance along x-axis)
y_adj <- 0.95   #location of regression equation (% distance along y-axis)
y_adj2 <- 0.05   #location of regression equation (% distance along y-axis)

#set flag variables for axis type for Fig3D and Fig3E
x_is_biovolume <- FALSE
y_is_biovolume <- FALSE
out_file <- "FF"

#open output plot file
if (pt_label){
   plot_labeled <- "(points ARE labeled).pdf"
} else {
   plot_labeled <- "(points NOT labeled).pdf"
}
plot_out <- paste0(pathFIG,"Figure 3.pdf")
pdf(file=plot_out, width=180/25.4, height=170/25.4)

#set global plotting options/parameters
op <- par(no.readonly = TRUE)
par(mfrow=c(3,2),
    mar=c(3,4,1.0,0.5),
    mgp=c(1.5,0.3,0),
    cex.axis=1.5,
    cex.lab=1.75,
    ps=6,
    family="sans")

#input volume
input_file <- paste0(pathDF,"Fig3F_rhodopsin_per_µm.csv")
Fig3F <- read.csv(file=input_file, header=TRUE)
Fig3F_subset <- unique(select(Fig3F, "genus","spherical.cell.volume..µm3."))

#input color code
input_file <- paste0(pathDF,"Genus_hex_color.csv")
gene_colors <- read.csv(file=input_file, header=TRUE)
gene_colors2 <- gene_colors[1:107,]
gene_colors2 <- rbind(gene_colors2,c("#000000","Unclassified","Unclassified"))
color_file <- paste0(pathDF,"Color_master_key.csv")
class_colors <- read.csv(file=color_file, header=TRUE)
gene_class_colors <- merge(x=gene_colors2, y=class_colors, 
                           by.x="class", by.y="Class", sort=FALSE, all=TRUE)

#input abundance data
input_file <- paste0(pathDF,"Fig2_abundance.csv")
input <- read.csv(file=input_file, header=TRUE)
#get unique genus and total abundance
unique_genes <- unique(input$genus)
unique_dates <- unique(input$date)
#reformat to a table with genes and then three days
data_table <- data.frame(unique_genes, 
                         as.numeric(input$cells_per_ml[input$date == unique_dates[1]]), 
                         as.numeric(input$cells_per_ml[input$date == unique_dates[2]]), 
                         as.numeric(input$cells_per_ml[input$date == unique_dates[3]]))
colnames(data_table) <- c("genus", unique_dates)
data_table2 <- merge(data_table, gene_class_colors, by.x="genus", by.y="GTDB_genus",
                      sort=FALSE, all=FALSE)

####################################################################################
# Figure 3A
####################################################################################
file_name <- paste0(pathDF,"Fig3_O2_genera_rRNA_diameter.csv")
Fig3_rRNA <- read.csv(file=file_name, header=TRUE)
allG <- merge(Fig3_rRNA, Fig3F_subset, by.x="genus", by.y="genus",
              sort=FALSE, all=FALSE)
allG <- merge(allG, gene_class_colors, by.x="genus", by.y="GTDB_genus",
              sort=FALSE, all=FALSE)

#keep only the unique entries (i.e. remove duplicates)
allG <- unique(allG)

#set plotting shapes
allG$shape <- rep(15, nrow(allG))
allG$shape[allG$date == "20190402"] <- 17
allG$shape[allG$date == "20190709"] <- 18

#keep only what's needed for plotting
plot_data <- select(allG, "Weighted_avg_O2_consumed_per_cell", "spherical.cell.volume..µm3.", 
                    "ClassColor", "shape", "genus")
colnames(plot_data) <- c("x", "y", "color", "shape", "genus")

#plotting parameters
x_lim <- c(1e-3,1e1)
x_main_label <- expression(paste("Respiration Rate (fmol O"[2]~"cell"^-1~"h"^-1*")"))
major_x_tic_location <- 10^c(-3:1)
major_x_tic_labels <- c("0.001", "0.01", "0.1", "1.0", "10")
minor_x_tic_location <- c(c(2:9)*1e-3, c(2:9)*1e-2, c(2:9)*1e-1, c(2:9)*1e0, c(2:9)*1e1)

y_lim <- c(0,0.4)
y_main_label <- expression(paste("Estimated Biovolume (",mu,"m"^3~"cell"^-1*")"))
major_y_tic_location <- seq(from=0, to=0.4, by=0.1)
major_y_tic_labels <- c("0.0 ", "0.1 ", "0.2 ", "0.3 ", "0.4 ")
minor_y_tic_location <- seq(from=0, to=0.4, by=0.2)

#set location of regression equation information
x_reg <- 10^(log10(x_lim[1]) + x_adj*(log10(x_lim[2])-log10(x_lim[1])))
y_reg <- y_lim[1] + y_adj*(y_lim[2]-y_lim[1])
reg_eqn_xy <- c(x_reg, y_reg)

label_points_xy <- rbind(c(0.003, 0.02, 0.01, 0.8),c(0.2, 0.01, 0.16, 0.12))
tab_y <- (0.4+0.075*(0.4+0))
tab_label <- "A"

par(mfg=c(1,1))
plot_logx(plot_data,
            x_lim, x_main_label, major_x_tic_location, major_x_tic_labels, minor_x_tic_location,
            y_lim, y_main_label, major_y_tic_location, major_y_tic_labels, minor_y_tic_location,
            reg_eqn_xy, pt_label, label_points_xy, tab_y, tab_label)


####################################################################################
# Figure 3C
####################################################################################
#keep only what's needed for plotting
plot_data <- select(allG, "spherical.cell.volume..µm3.", "X16S_copies_per_cell", 
                    "ClassColor", "shape", "genus")
colnames(plot_data) <- c("x", "y", "color", "shape", "genus")

#plotting parameters
x_lim <- c(0,0.4)
x_main_label <- expression(paste("Estimated Biovolume (",mu,"m"^3~"cell"^-1*")"))
major_x_tic_location <- c(0.0, 0.1, 0.2, 0.3, 0.4)
major_x_tic_labels <- c("0.0", "0.1", "0.2", "0.3", "0.4")
minor_x_tic_location <- seq(from=0, to=0.4, by=0.02)

y_lim <- c(1e1,1e4)
#units <- expression(paste("(cell"^-1*")"))
#mtext(text=units, side=2, line=1.5)

y_main_label <- expression(paste("Count of 16S rRNA Transcripts (cell"^-1*")"))
major_y_tic_location <- 10^(1:4)
major_y_tic_labels <- c("10 ", "100 ", "1000 ", "10000 ")
minor_y_tic_location <- c(c(2:9)*1e1, c(2:9)*1e2, c(2:9)*1e3, c(2:9)*1e4)

#set location of regression equation information
x_reg <- x_lim[1] + x_adj*(x_lim[2]-x_lim[1])
y_reg <- 10^(log10(y_lim[1]) + y_adj2*(log10(y_lim[2])-log10(y_lim[1])))
reg_eqn_xy <- c(x_reg, y_reg)

#set location of point labels
label_points_xy <- rbind(c(0.08, 0.15, 0.04, 0.14),c(40, 100, 20, 5000))
tab_y <- 10^(4+0.075*(4-1))
tab_label <- "C"

par(mfg=c(2,1), ps=6)
plot_logy(plot_data,
            x_lim, x_main_label, major_x_tic_location, major_x_tic_labels, minor_x_tic_location,
            y_lim, y_main_label, major_y_tic_location, major_y_tic_labels, minor_y_tic_location,
            reg_eqn_xy, pt_label, label_points_xy, tab_y, tab_label)


####################################################################################
# Figure 3D
####################################################################################
#input data and create plotting data
file_name <- paste0(pathDF,"Fig3E_rRNA_O2_rhodhopsin.csv")
Fig3E_rRNA <- read.csv(file=file_name, header=TRUE)
Fig3E_rRNA <- unique(select(Fig3E_rRNA, "genus", "date", 
                     "X16S_copies_per_cell",
                     "Weighted_avg_O2_consumed_per_cell", 
                     "proteorhodopsin_transcripts_per_cell"))
allG <- merge(Fig3E_rRNA, Fig3F_subset, by.x="genus", by.y="genus",
              sort=FALSE, all=FALSE)
allG <- merge(allG, gene_class_colors, by.x="genus", by.y="GTDB_genus",
              sort=FALSE, all=FALSE)

#keep only the unique entries (i.e. remove duplicates)
allG <- unique(allG)

#keep only data within desired range
keep <- (allG$proteorhodopsin_transcripts_per_cell > min_proteorhodopsin) & 
   (allG$proteorhodopsin_transcripts_per_cell < max_proteorhodopsin)
allG <- allG[keep,]

#calculate quantites based on biovolume
allG$count_per_volume <- allG$X16S_copies_per_cell/allG$spherical.cell.volume..µm3.
allG$O2_per_volume <- allG$Weighted_avg_O2_consumed_per_cell/allG$spherical.cell.volume..µm3.

#set plotting shapes
allG$shape <- rep(15, nrow(allG))
allG$shape[allG$date == "20190402"] <- 17
allG$shape[allG$date == "20190709"] <- 18

if (x_is_biovolume & y_is_biovolume){
   #keep only what's needed for plotting
   plot_data <- select(allG, "O2_per_volume", "count_per_volume", "ClassColor", "shape", "genus")
   colnames(plot_data) <- c("x", "y", "color", "shape", "genus")
   
   #plotting parameters
   x_lim <- c(0.001,20)
   x_main_label <- expression(paste("Respiration Rate (fmol O"[2]," ",mu,"m"^-3," h"^-1*")"))
   major_x_tic_location <- c(1e-3,1e-2,1e-1,1e0,1e1)
   major_x_tic_labels <- c("0.001","0.01", "0.1", "1.0", "10")
   minor_x_tic_location <- c(c(2:9)*1e-3,c(2:9)*1e-2, c(2:9)*1e-1, c(2:9)*1e0,20)
   
   y_lim <- c(1,1e5)
   y_main_label <- expression(paste("Count of 16S rRNA Transcripts (",mu,"m"^-3*")"))
   major_y_tic_location <- 10^(0:5)
   major_y_tic_labels <- c("1 ","10 ", "100 ", "1000 ", "10000 ", "100000 ")
   minor_y_tic_location <- c(c(2:9)*1e0, c(2:9)*1e1, c(2:9)*1e2, c(2:9)*1e3, c(2:9)*1e4)
   
   #set location of regression equation information
   x_reg <- 10^(log10(x_lim[1]) + x_adj*(log10(x_lim[2])-log10(x_lim[1])))
   y_reg <- 10^(log10(y_lim[1]) + y_adj*(log10(y_lim[2])-log10(y_lim[1])))
   reg_eqn_xy <- c(x_reg, y_reg)
   
   #set location of point labels
   label_points_xy <- rbind(c(1, 0.2, 0.08, 2),c(1000,200,60,90000))
   tab_y <- 10^(5+0.075*(6))
} else if(x_is_biovolume & !y_is_biovolume){
   #keep only what's needed for plotting
   plot_data <- select(allG, "O2_per_volume", "X16S_copies_per_cell", 
                       "ClassColor", "shape", "genus")
   colnames(plot_data) <- c("x", "y", "color", "shape", "genus")
   
   #plotting parameters
   x_lim <- c(0.001,20)
   x_main_label <- expression(paste("Respiration Rate (fmol O"[2]," ",mu,"m"^-3," h"^-1*")"))
   major_x_tic_location <- c(1e-3,1e-2,1e-1,1e0,1e1)
   major_x_tic_labels <- c("0.001","0.01", "0.1", "1.0", "10")
   minor_x_tic_location <- c(c(2:9)*1e-3,c(2:9)*1e-2, c(2:9)*1e-1, c(2:9)*1e0,20)
   
   y_lim <- c(1,10000)
   y_main_label <- expression(paste("Count of 16S rRNA Transcripts (cell"^-1*")"))
   major_y_tic_location <- 10^(0:4)
   major_y_tic_labels <- c("1 ", "10 ","100 ", "1000 ", "10000 ")
   minor_y_tic_location <- c(c(2:9)*1e0, c(2:9)*1e1, c(2:9)*1e2, c(2:9)*1e3)
   
   #set location of regression equation information
   x_reg <- 10^(log10(x_lim[1]) + x_adj*(log10(x_lim[2])-log10(x_lim[1])))
   y_reg <- 10^(log10(y_lim[1]) + y_adj*(log10(y_lim[2])-log10(y_lim[1])))
   reg_eqn_xy <- c(x_reg, y_reg)
   
   #set location of point labels
   label_points_xy <- rbind(c(0.1, 2, 0.05, 1.95),c(20,100,60,800))
   tab_y <- 10^(4+0.075*(3))
} else if(!x_is_biovolume & y_is_biovolume){
   #keep only what's needed for plotting
   plot_data <- select(allG, "Weighted_avg_O2_consumed_per_cell", "count_per_volume", 
                       "ClassColor", "shape", "genus")
   colnames(plot_data) <- c("x", "y", "color", "shape", "genus")
   
   #plotting parameters
   x_lim <- c(0.001,10)
   x_main_label <- expression(paste("Respiration Rate (fmol O"[2]," cell"^-1," h"^-1*")"))
   major_x_tic_location <- c(1e-3,1e-2,1e-1,1e0,1e1)
   major_x_tic_labels <- c("0.001","0.01", "0.1", "1.0", "10")
   minor_x_tic_location <- c(c(2:9)*1e-3,c(2:9)*1e-2, c(2:9)*1e-1, c(2:9)*1e0)
   
   y_lim <- c(1,1e5)
   y_main_label <- expression(paste("Count of 16S rRNA Transcripts (",mu,"m"^-3*")"))
   major_y_tic_location <- 10^(0:5)
   major_y_tic_labels <- c("1 ","10 ", "100 ", "1000 ", "10000 ", "100000 ")
   minor_y_tic_location <- c(c(2:9)*1e0, c(2:9)*1e1, c(2:9)*1e2, c(2:9)*1e3, c(2:9)*1e4)
   
   #set location of regression equation information
   x_reg <- 10^(log10(x_lim[1]) + x_adj*(log10(x_lim[2])-log10(x_lim[1])))
   y_reg <- 10^(log10(y_lim[1]) + y_adj*(log10(y_lim[2])-log10(y_lim[1])))
   reg_eqn_xy <- c(x_reg, y_reg)
   
   #set location of point labels
   label_points_xy <- rbind(c(0.01, 0.01, 0.05, .4),c(10,40000,60,600))
   tab_y <- 10^(5+0.075*(5))
} else if(!x_is_biovolume & !y_is_biovolume){
   #keep only what's needed for plotting
   plot_data <- select(allG, "Weighted_avg_O2_consumed_per_cell", "X16S_copies_per_cell", 
                       "ClassColor", "shape", "genus")
   colnames(plot_data) <- c("x", "y", "color", "shape", "genus")
   
   #plotting parameters
   x_lim <- c(0.001,10)
   x_main_label <- expression(paste("Respiration Rate (fmol O"[2]," cell"^-1," h"^-1*")"))
   major_x_tic_location <- c(1e-3,1e-2,1e-1,1e0,1e1)
   major_x_tic_labels <- c("0.001","0.01", "0.1", "1.0", "10")
   minor_x_tic_location <- c(c(2:9)*1e-3,c(2:9)*1e-2, c(2:9)*1e-1, c(2:9)*1e0)
   
   y_lim <- c(1,10000)
   y_main_label <- expression(paste("Count of 16S rRNA Transcripts (cell"^-1*")"))
   major_y_tic_location <- 10^(0:4)
   major_y_tic_labels <- c("1 ", "10 ","100 ", "1000 ", "10000 ")
   minor_y_tic_location <- c(c(2:9)*1e0, c(2:9)*1e1, c(2:9)*1e2, c(2:9)*1e3)
   
   #set location of regression equation information
   x_reg <- 10^(log10(x_lim[1]) + x_adj*(log10(x_lim[2])-log10(x_lim[1])))
   y_reg <- 10^(log10(y_lim[1]) + y_adj*(log10(y_lim[2])-log10(y_lim[1])))
   reg_eqn_xy <- c(x_reg, y_reg)
   
   #set location of point labels
   label_points_xy <- rbind(c(0.01, 0.004, 0.05, .4),c(10,4000,60,600))
   tab_y <- 10^(4+0.075*(3))
} 

tab_label <- "D"
par(mfg=c(2,2), ps=6)
plot_logxy(plot_data,
           x_lim, x_main_label, major_x_tic_location, major_x_tic_labels, minor_x_tic_location,
           y_lim, y_main_label, major_y_tic_location, major_y_tic_labels, minor_y_tic_location,
           reg_eqn_xy, pt_label, label_points_xy, tab_y, tab_label)

####################################################################################
# Figure 3E
####################################################################################
file_name <- paste0(pathDF,"Fig3_O2_genera_mRNA_diameter.csv")
Fig3E_mRNA <- read.csv(file=file_name, header=TRUE)

allG <- merge(Fig3E_mRNA, Fig3F_subset, by.x="genus", by.y="genus",
              sort=FALSE, all=FALSE)
allG <- merge(allG, gene_class_colors, by.x="genus", by.y="GTDB_genus",
              sort=FALSE, all=FALSE)

#keep only the unique entries (i.e. remove duplicates)
allG <- unique(allG)

#calculate quantites based on biovolume
allG$count_per_volume <- allG$mRNA_transcripts_per_cell/allG$spherical.cell.volume..µm3.
allG$O2_per_volume <- allG$Weighted_avg_O2_consumed_per_cell/allG$spherical.cell.volume..µm3.

#set plotting shapes
allG$shape <- rep(15, nrow(allG))
allG$shape[allG$date == "20190402"] <- 17
allG$shape[allG$date == "20190709"] <- 18

if (x_is_biovolume & y_is_biovolume){
   #keep only what's needed for plotting
   plot_data <- select(allG, "O2_per_volume", "count_per_volume", "ClassColor", "shape", "genus")
   colnames(plot_data) <- c("x", "y", "color", "shape", "genus")
   
   #plotting parameters
   x_lim <- c(1e-2,2e1)
   x_main_label <- expression(paste("Respiration Rate (fmol O"[2]," ",mu,"m"^-3," h"^-1*")"))
   major_x_tic_location <- c(0.01, 0.1, 1.0, 10.)
   major_x_tic_labels <- c("0.01", "0.1", "1.0", "10")
   minor_x_tic_location <- c(c(2:9)*1e-2, c(2:9)*1e-1, c(2:9)*1e0,20)
   
   y_lim <- c(1e2,1e4)
   y_main_label <- expression(paste("Count of 16S mRNA Transcripts (",mu,"m"^-3*")"))
   major_y_tic_location <- 10^(2:4)
   major_y_tic_labels <- c("100 ","1000 ","10000 ")
   minor_y_tic_location <- c(c(2:9)*1e2, c(2:9)*1e3)
   
   #set location of regression equation information
   x_reg <- 10^(log10(x_lim[1]) + x_adj*(log10(x_lim[2])-log10(x_lim[1])))
   y_reg <- 10^(log10(y_lim[1]) + y_adj*(log10(y_lim[2])-log10(y_lim[1])))
   reg_eqn_xy <- c(x_reg, y_reg)
   
   #set location of point labels
   label_points_xy <- rbind(c(0.5, 0.08, 0.08, 2),c(2250,6000,130,500))
   tab_y <- 10^(4+0.075*(2))
} else if(x_is_biovolume & !y_is_biovolume){
   #keep only what's needed for plotting
   plot_data <- select(allG, "O2_per_volume", "mRNA_transcripts_per_cell", "ClassColor", "shape", "genus")
   colnames(plot_data) <- c("x", "y", "color", "shape", "genus")
   
   #plotting parameters
   x_lim <- c(1e-2,2e1)
   x_main_label <- expression(paste("Respiration Rate (fmol O"[2]," ",mu,"m"^-3," h"^-1*")"))
   major_x_tic_location <- c(0.01, 0.1, 1.0, 10.)
   major_x_tic_labels <- c("0.01", "0.1", "1.0", "10")
   minor_x_tic_location <- c(c(2:9)*1e-2, c(2:9)*1e-1, c(2:9)*1e0,20)
   
   y_lim <- c(1,1000)
   y_main_label <- expression(paste("Count of 16S mRNA Transcripts (cell"^-1*")"))
   major_y_tic_location <- c(1,10,100,1000)
   major_y_tic_labels <- c("1 ", "10 ", "100 ","1000 ")
   minor_y_tic_location <- c(c(2:9)*1e0, c(2:9)*1e1, c(2:9)*1e2)
   
   #set location of regression equation information
   x_reg <- 10^(log10(x_lim[1]) + x_adj*(log10(x_lim[2])-log10(x_lim[1])))
   y_reg <- 10^(log10(y_lim[1]) + y_adj*(log10(y_lim[2])-log10(y_lim[1])))
   reg_eqn_xy <- c(x_reg, y_reg)
   
   #set location of point labels
   label_points_xy <- rbind(c(0.1, 0.1, 0.2, 2),c(5,600,3,5))
   tab_y <- 10^(3+0.075*(3))
} else if(!x_is_biovolume & y_is_biovolume){
   #keep only what's needed for plotting
   plot_data <- select(allG, "Weighted_avg_O2_consumed_per_cell", "count_per_volume", 
                       "ClassColor", "shape", "genus")
   colnames(plot_data) <- c("x", "y", "color", "shape", "genus")
   
   #plotting parameters
   x_lim <- c(1e-3,1e1)
   x_main_label <- expression(paste("Respiration Rate (fmol O"[2]," cell"^-1," h"^-1*")"))
   major_x_tic_location <- c(0.001, 0.01, 0.1, 1.0, 10.)
   major_x_tic_labels <- c("0.001", "0.01", "0.1", "1.0", "10")
   minor_x_tic_location <- c(c(2:9)*1e-3, c(2:9)*1e-2, c(2:9)*1e-1, c(2:9)*1e0)
   
   y_lim <- c(1e2,1e4)
   y_main_label <- expression(paste("Count of 16S mRNA Transcripts (",mu,"m"^-3*")"))
   major_y_tic_location <- 10^(2:4)
   major_y_tic_labels <- c("100 ","1000 ","10000 ")
   minor_y_tic_location <- c(c(2:9)*1e2, c(2:9)*1e3)
   
   #set location of regression equation information
   x_reg <- 10^(log10(x_lim[1]) + x_adj*(log10(x_lim[2])-log10(x_lim[1])))
   y_reg <- 10^(log10(y_lim[1]) + y_adj*(log10(y_lim[2])-log10(y_lim[1])))
   reg_eqn_xy <- c(x_reg, y_reg)
   
   #set location of point labels
   label_points_xy <- rbind(c(0.01, 0.003, 0.02, 0.03),c(2000,6000,1400,180))
   tab_y <- 10^(5+0.075*(4))
} else if(!x_is_biovolume & !y_is_biovolume){
   #keep only what's needed for plotting
   plot_data <- select(allG, "Weighted_avg_O2_consumed_per_cell", "mRNA_transcripts_per_cell", 
                       "ClassColor", "shape", "genus")
   colnames(plot_data) <- c("x", "y", "color", "shape", "genus")
   
   #plotting parameters
   x_lim <- c(1e-3,1e1)
   x_main_label <- expression(paste("Respiration Rate (fmol O"[2]," cell"^-1," h"^-1*")"))
   major_x_tic_location <- c(0.001, 0.01, 0.1, 1.0, 10.)
   major_x_tic_labels <- c("0.001", "0.01", "0.1", "1.0", "10")
   minor_x_tic_location <- c(c(2:9)*1e-3, c(2:9)*1e-2, c(2:9)*1e-1, c(2:9)*1e0)
   
   y_lim <- c(1,1000)
   y_main_label <- expression(paste("Count of mRNA Transcripts (cell"^-1*")"))
   major_y_tic_location <- c(1,10,100,1000)
   major_y_tic_labels <- c("1 ", "10 ", "100 ","1000 ")
   minor_y_tic_location <- c(c(2:9)*1e0, c(2:9)*1e1, c(2:9)*1e2)
   
   #set location of regression equation information
   x_reg <- 10^(log10(x_lim[1]) + x_adj*(log10(x_lim[2])-log10(x_lim[1])))
   y_reg <- 10^(log10(y_lim[1]) + y_adj*(log10(y_lim[2])-log10(y_lim[1])))
   reg_eqn_xy <- c(x_reg, y_reg)
   
   #set location of point labels
   label_points_xy <- rbind(c(0.01, 0.0025, 0.0035, 0.35),c(6,600,3,10))
   tab_y <- 10^(3+0.075*(3))
}

tab_label <- "E"
par(mfg=c(3,1), ps=6)
plot_logxy(plot_data,
           x_lim, x_main_label, major_x_tic_location, major_x_tic_labels, minor_x_tic_location,
           y_lim, y_main_label, major_y_tic_location, major_y_tic_labels, minor_y_tic_location,
           reg_eqn_xy, pt_label, label_points_xy, tab_y, tab_label)

####################################################################################
# Figure 3B
####################################################################################
input_file <- paste0(pathDF,"All_days_avg_O2_consumption_and_abundance_virus_percentage.csv")
percent_data <- read.csv(file=input_file, header=TRUE)
keep <- !is.na(percent_data$Weighted_avg_O2_consumed_per_cell)
percent_data <- percent_data[keep,]
allG <- merge(percent_data, Fig3F_subset, by.x="genus", by.y="genus",
              sort=FALSE, all=FALSE)
allG <- merge(allG, gene_class_colors, by.x="genus", by.y="GTDB_genus",
              sort=FALSE, all=FALSE)

#Growth rate plot
input_file <- paste0(pathDF,"growth_rate.csv")
growth_rate <- read.csv(file=input_file, header=TRUE)
growth_rate <- growth_rate[!is.na(growth_rate$Average_doubling_time),]
abc <- unique(merge(x=growth_rate, y=allG, 
             by.x="GTDB_genus", by.y="genus", 
             all=FALSE, sort=FALSE))
abc <- abc[abc$Weighted_avg_O2_consumed_per_cell > min_O2_keep,]
abc$shape <- rep(15, nrow(abc))
abc$shape[abc$date == "20190402"] <- 17
abc$shape[abc$date == "20190709"] <- 18
plot_data <- select(abc, "Weighted_avg_O2_consumed_per_cell", "Average_doubling_time",
                    "ClassColor", "shape", "GTDB_genus")
colnames(plot_data) <- c("x", "y", "color", "shape", "genus")

#plotting parameters
x_lim <- c(1e-3,1e1)
x_main_label <- expression(paste("Respiration Rate (fmol O"[2]~"cell"^-1~"h"^-1*")"))
major_x_tic_location <- 10^c(-3:1)
major_x_tic_labels <- c("0.001", "0.01", "0.1", "1.0", "10")
minor_x_tic_location <- c(c(2:9)*1e-3, c(2:9)*1e-2, c(2:9)*1e-1, c(2:9)*1e0, c(2:9)*1e1)

y_lim <- c(0,12)
y_main_label <- "Predicted Minimum Doubling Time (h)"
major_y_tic_location <- seq(from=0, to=12, by=2)
major_y_tic_labels <- c("0 ", "2 ", "4 ", "6 ", "8 ", "10 ", "12 ")
minor_y_tic_location <- seq(from=1, to=11, by=2)

#set location of regression equation information
x_reg <- 10^(log10(x_lim[1]) + x_adj*(log10(x_lim[2])-log10(x_lim[1])))
y_reg <- y_lim[1] + y_adj*(y_lim[2]-y_lim[1])
reg_eqn_xy <- c(x_reg, y_reg)

#set location of point labels
label_points_xy <- rbind(c(0.01, 0.03, 0.1, 0.4),c(9, 7.5, 5.5, 4.5))
tab_y <- (12+0.075*(12+0))
tab_label <- "B"

par(mfg=c(1,2), ps=6)
plot_logx(plot_data,
            x_lim, x_main_label, major_x_tic_location, major_x_tic_labels, minor_x_tic_location,
            y_lim, y_main_label, major_y_tic_location, major_y_tic_labels, minor_y_tic_location,
            reg_eqn_xy, pt_label, label_points_xy, tab_y, tab_label)
####################################################################################
# Figure LEGEND
####################################################################################
#add the legend in the bottom right space
#first, draw an empty graph
par(mfg=c(3,2), ps=6)
plot(x=NA, y=NA, xlim=c(0,1), ylim=c(0,1), 
     xlab="", ylab="",bty="n", xaxt="n", yaxt="n")

#add legend for all plots at bottom right of page
legend_info <- unique(gene_class_colors[,c(1,4)])
rownames(legend_info) <- 1:nrow(legend_info)
legend(x=0.0, y=1.1, ncol=1, pch=22, 
       legend=legend_info$class[c(6,9)],
       col=legend_info$ClassColor[c(6,9)], 
       pt.bg=legend_info$ClassColor[c(6,9)],
       xpd=NA, bty="n", y.intersp = 0.7, 
       title="Archaeal Classes", title.adj=c(0,0), cex=1.75)
legend(x=0.0, y=0.7, ncol=1, pch=22, 
       legend=legend_info$class[c(1,3,5,4,10,11,2)],
       col=legend_info$ClassColor[c(1,3,5,4,10,11,2)], 
       pt.bg=legend_info$ClassColor[c(1,3,5,4,10,11,2)],
       xpd=NA, bty="n", y.intersp = 0.7, 
       title="Bacterial Classes", title.adj=c(0,0), cex=1.75)
legend(x=0.5, y=1.1, legend=rbind("Oct 30, 2018","Apr 02, 2019","Jul 09, 2019"),
       pch=c(22,23,24), col="black", pt.bg="black",
       xpd=NA, bty="n", y.intersp = 0.7, 
       title="Sampling Date", title.adj=c(0,0), cex=1.75, pt.cex=1.5)

par(op)
dev.off()