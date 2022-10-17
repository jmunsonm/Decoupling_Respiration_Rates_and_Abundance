rm(list = ls())

library(plot3D)
library(dplyr)
library(viridis)

path0 <- "~/Desktop/Git_upload/"
pathDF <- "~/Desktop/Git_upload/Data Files/"
pathFIG <- "~/Desktop/Git_upload/Figures/"

delta_y <- 0.02
ylim1 <- 0.10
ylim5 <- 0.98
y_height <- (ylim5-ylim1)/4-2*delta_y
ylim2 <- ylim1+y_height+delta_y
ylim3 <- ylim2+y_height+delta_y
ylim4 <- ylim3+y_height+3*delta_y

input_file <- paste0(pathDF,"data_for_fig2B.csv")
input <- read.csv(file=input_file, header=TRUE)

#change divergence into numerical values for plot color
unique_div <- c("Different Class", "Same Class", "Same Order", "Same Family", "Same Genus")
plot_col <- rep(NA, nrow(input))
for (i in seq_along(unique_div)){
  change <- (input$divergence == unique_div[i])
  plot_col[change] <- i
}

#data for the line of median values
x_unique <- unique(input$AAI_round)
x_unique <- x_unique[order(x_unique, decreasing=FALSE)]
x_median <- matrix(NA, ncol=3, nrow=length(x_unique))
for (i in seq_along(x_unique)){
  x_rows <- (input$AAI_round == x_unique[i])
  x_data <- input$O2_ratio[x_rows]
  x_quantile <- quantile(x_data, c(.10,.50,.90))
  x_median[i,] <- x_quantile
}
x_poly <- c(x_unique, x_unique[length(x_unique):1])
y_poly <- c(x_median[,1], x_median[(length(x_unique):1),3])

pdf(file=paste0(pathFIG,"Figure 2.pdf"), width=90/25.4, height=247/25.4)
op <- par(no.readonly = TRUE)
par(mar=c(0,0,0,0), cex = 1, cex.main = 1)
par(ps=6, fig=c(0.15,0.99,ylim4+delta_y,ylim5))

#basic plot
points2D(x=input$Mean.AAI, y=input$O2_ratio, 
         ylim=c(1,2e3), yaxt='n', ylab="",
         xlim=c(30,100), xlab="", xaxt='n', 
         pch=16, cex=0.5,
         colvar=plot_col, col=viridis(5),
         colkey=FALSE,
         log="y")
polygon2D(x=x_poly, y=y_poly, 
          yaxt='n', log="y",
          col=rgb(10, 10, 10, max = 255, alpha = 140),
          add=TRUE)
lines2D(x=x_unique, y=x_median[,2], col="black", 
        yaxt='n', log="y", lwd=2,
        add=TRUE)
legend(x="topright",legend=unique_div, pch=rep(16, length(unique_div)), 
       col=viridis(5),
       cex=1.0, y.intersp = 0.5)
#x-axis annotation
mtext(text="Pairwise AAI", side=1, line=0.5, xpd=NA)
axis(side=1,
     at=seq(from=30, to=100, by=10),labels=FALSE)
axis(labels=seq(from=30, to=100, by=10), 
     at=seq(from=30, to=100, by=10),
     side=1, xpd=NA, tick=FALSE, padj=-3)
#y-axis annotation
mtext(text=expression("O"[2]*" Consumption Ratio"),
      side=2, line=1.75, xpd=NA, srt=90)
axis(side=2, 
     at=c(1e0, 1e1, 1e2, 1e3, 1e4), 
     labels=c(expression("1x10"^0), expression("1x10"^1), 
              expression("1x10"^2), expression("1x10"^3), expression("1x10"^4)),
     las=1, hadj=0.6, padj=0.40)
axis(side=2, 
     at=c(c(2:9)*1e0, c(2:9)*1e1, c(2:9)*1e2),  #, c(2:9)*1e3), 
     labels=FALSE, tcl=-0.25)

par(ps=8)
text(x=14.5, y=3e3, labels="B", cex=1, adj=c(0,0), xpd=NA, font=2)
par(ps=6)

#####################################################################################

#input color code
input_file <- paste0(pathDF,"Genus_hex_color.csv")
gene_colors <- read.csv(file=input_file, header=TRUE)
gene_colors <- rbind(gene_colors[1:107,],c("#000000","Unclassified","Unclassified"))

combined_names <- c("Pelagibacter", 
                    "SCGC-AAA076-P13", 
                    "MS024-2A", 
                    "D2472", 
                    "Thioglobus", 
                    "Hel1-33-131",    
                    "BACL14", 
                    "Planktomarina", 
                    "Amylibacter", 
                    "MAG-121220-bin8", 
                    "Polaribacter", 
                    "SCGC-AAA160-P02",
                    "HC6-5", 
                    "GCA-002733185", 
                    "UBA7428", 
                    "ASP10-02a", 
                    "UBA4465", 
                    "LFER01", 
                    "Unclassified")
v_color <- c("red","blue1","darkgoldenrod4")
other_color <- "black"
name_color <- c(rep(v_color[1],7),rep(v_color[2],2),v_color[1],rep(v_color[3],8),"black")
combined_names <- cbind(combined_names, name_color)
colnames(combined_names) <- c("genus", "name_color")

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
colnames(data_table) <- c("genus",unique_dates)
#add total abundance
total_abund <- data_table[,2] + data_table[,3] + data_table[,4]
data_table2 <- data.frame(data_table,total_abund)
data_table2b <- merge(data_table2, gene_colors[1:108,], by.x="genus", by.y="GTDB_genus",
                      sort=FALSE, all=FALSE)

#add colors by class
input_file <- paste0(pathDF,"Color_master_key.csv")
class_color <- read.csv(input_file, header=TRUE)
class_color <- rbind(class_color,c("Unclassified","#000000"))
data_table3 <- merge(x=data_table2b, y=class_color, by.x="class", by.y="Class",
                     sort=FALSE, all=FALSE)#sort by total abundance
new_order <- order(data_table3$total_abund, decreasing=TRUE)
data_table3 <- data_table3[new_order,]

#first plot
#basic plot
par(fig=c(0.15,0.99,ylim3,ylim3+y_height), new=TRUE)
par(ps=6)
y_data <- merge(x=combined_names, y=data_table3, by.x="genus", by.y="genus",
                sort=FALSE, all=FALSE)
other_1 <- 1-sum(y_data$X20181030, rm.na=TRUE)/sum(data_table3$X20181030, rm.na=TRUE)
other_2 <- 1-sum(y_data$X20190402, rm.na=TRUE)/sum(data_table3$X20190402, rm.na=TRUE)
other_3 <- 1-sum(y_data$X20190709, rm.na=TRUE)/sum(data_table3$X20190709, rm.na=TRUE)

unclass_1 <- y_data$X20181030[19]/sum(data_table3$X20181030, rm.na=TRUE)
unclass_2 <- y_data$X20190402[19]/sum(data_table3$X20190402, rm.na=TRUE)
unclass_3 <- y_data$X20190709[19]/sum(data_table3$X20190709, rm.na=TRUE)

other_abund <- c(other_1+unclass_1, other_2+unclass_2, other_3+unclass_3)

x_seq <- seq(from=2, to=74, by=4)
points2D(x=x_seq[1:18]-1, 
         y=(100*y_data$X20181030/sum(data_table3$X20181030, rm.na=TRUE))[1:18], 
         ylim=c(1e-2,1e2), 
         xlim=c(1,4*nrow(y_data)-2),
         xaxt='n', yaxt='n', pch=15, xlab="", ylab="", 
         col=y_data$ClassColor, log="y", cex=0.4)
points2D(x=x_seq[1:18], 
         y=(100*y_data$X20190402/sum(data_table3$X20190402, rm.na=TRUE))[1:18], 
         xaxt='n', yaxt='n', pch=15, 
         col=y_data$ClassColor, log="y", cex=0.4, add=TRUE)
points2D(x=x_seq[1:18]+1, 
         y=(100*y_data$X20190709/sum(data_table3$X20190709, rm.na=TRUE))[1:18], 
         xaxt='n', yaxt='n', pch=15, 
         col=y_data$ClassColor, log="y", cex=0.4, add=TRUE)
points2D(x=c(4*nrow(y_data)-c(3,2,1)), 
         y=100*other_abund, 
         xaxt='n', yaxt='n', pch=15, 
         col=other_color, log="y", cex=0.4, add=TRUE)
#add grey vertical lines between genus
for (iline in seq(from=0, to=79, by=4)){
  lines(x=c(iline,iline), y=c(0.01,100), type="l", col="grey50", lwd=0.5) 
}
#x-axis annotation
#text(x=x_seq, y=0.32e1, labels=y_data$genus, 
#     srt=90, adj=c(1,0.5), cex=0.6, xpd=NA, col=y_data$name_color)
axis(side=1, at=c(x_seq,78), labels=FALSE)
#y-axis annotation
mtext(text="Cells per ml (daily %)", 
      side=2, line=1.75, xpd=NA, srt=90)
axis(side=2, 
     at=c(1e-2,1e-1,1e0,1e1,1e2), 
     labels=c(expression("1x10"^-2), expression("1x10"^-1), 
              expression("1x10"^0), expression("1x10"^1), expression("1x10"^2)),
     las=1, hadj=0.6, padj=0.40)
axis(side=2, 
     at=c(c(2:9)*1e-2, c(2:9)*1e-1, c(2:9)*1e0, c(2:9)*1e1), 
     labels=FALSE, tcl=-0.25)

#top x-axis annotation (date)
axis(side=3, at=c(1,3), labels=FALSE, tck=-0.02)
axis(side=3, at=c(2), labels=FALSE, tck=-0.07)
text(x=1, y=2.3e2, labels="30 Oct 18", adj=c(0.95,0.5), xpd=NA, font=1)
text(x=2, y=3.5e2, labels="02 Apr 19", adj=c(0.5,0), xpd=NA, font=1)
text(x=3, y=2.3e2, labels="09 Jul 19", adj=c(0,0.5), xpd=NA, font=1)

par(ps=8)
text(x=-16., y=2e2, labels="C", cex=1, adj=c(0,0), xpd=NA, font=2)

#second plot
input_file <- paste0(pathDF,"per_cell_O2_consumption_rates_for_fiveSAG_genera.csv")
input <- read.csv(file=input_file, header=TRUE)
box_data2 <- matrix(NA, nrow=500, ncol=4*nrow(y_data)+3)
#add respiration rate to correct genus
unused_resp <- input$Predicted_respiration_rate_fmol_O2perhr
for (i in 1:nrow(y_data)){
  col1 <- 4*(i-1) + 1
  col2 <- col1 + 1
  col3 <- col2 + 1
  keep <- (y_data$genus[i] == input$GTDB_classification) & (input$FACS_date == 181030)
  unused_resp[keep] <- NA
  if(sum(keep, na.rm=TRUE) >= 1){
    resp_values <- input$Predicted_respiration_rate_fmol_O2perhr[keep]
    if(length(resp_values) >= 1){
      box_data2[1:length(resp_values),col1] <- resp_values
    }
  }
  keep <- (y_data$genus[i] == input$GTDB_classification) & (input$FACS_date == 190402)
  unused_resp[keep] <- NA
  if(sum(keep, na.rm=TRUE) >= 1){
    resp_values <- input$Predicted_respiration_rate_fmol_O2perhr[keep]
    if(length(resp_values) >= 1){
      box_data2[1:length(resp_values),col2] <- resp_values
    }
  }
  
  keep <- (y_data$genus[i] == input$GTDB_classification) & (input$FACS_date == 190709)
  unused_resp[keep] <- NA
  if(sum(keep, na.rm=TRUE) >= 1){
    resp_values <- input$Predicted_respiration_rate_fmol_O2perhr[keep]
    if(length(resp_values) >= 1){
      box_data2[1:length(resp_values),col3] <- resp_values
    }
  }
}

unused1 <- !is.na(unused_resp) & (input$FACS_date == 181030)
box_data2[1:sum(unused1),(1+(4*nrow(y_data)))] <- unused_resp[unused1]

unused2 <- !is.na(unused_resp) & (input$FACS_date == 190402)
box_data2[1:sum(unused2),(2+(4*nrow(y_data)))] <- unused_resp[unused2]

unused3 <- !is.na(unused_resp) & (input$FACS_date == 190709)
box_data2[1:sum(unused3),(3+(4*nrow(y_data)))] <- unused_resp[unused3]

#combine unclassified and other into one
used_73 <- sum(!is.na(box_data2[,73]))
used_74 <- sum(!is.na(box_data2[,74]))
used_75 <- sum(!is.na(box_data2[,75]))
used_77 <- sum(!is.na(box_data2[,77]))
used_78 <- sum(!is.na(box_data2[,78]))
used_79 <- sum(!is.na(box_data2[,79]))
box_data2[1:(used_73+used_77),73] <- c(box_data2[1:used_73,73],box_data2[1:used_77,77])
box_data2[1:(used_74+used_78),74] <- c(box_data2[1:used_74,74],box_data2[1:used_78,78])
box_data2[1:(used_75+used_79),75] <- c(box_data2[1:used_75,75],box_data2[1:used_79,79])

#basic plot
par(fig=c(0.15,0.99,ylim2,ylim2+y_height), new=TRUE)
par(ps=6)
color1 <- c(y_data$ClassColor,other_color)
color2 <- rbind(color1,color1,color1,"#FFFFFF")
box_col <- matrix(color2, ncol=1, byrow=FALSE)
boxplot(x=box_data2[,1:75], #names=c(y_data$genus,"All Others"), 
        col=box_col,   #c(y_data$ClassColor,other_color), 
        log="y", ylim=c(0.001,100),  
        xlim=c(1,4*nrow(y_data)-2),
        xaxt='n', yaxt="n",
        pars=list(boxwex=0.8, staplewex=0.2, outwex=0.2, 
                  boxlwd=0.2, medlwd=0.2, whisklwd=0.2, outlwd=0.2,
                  outcex=0.2, outbg=box_col, outpch=19,
                  outcol=box_col))
#add grey vertical lines between genus
for (iline in seq(from=0, to=79, by=4)){
  lines(x=c(iline,iline), y=c(0.001,100), type="l", col="grey50", lwd=0.5)
}
#x-axis annotation
axis(side=1, at=seq(from=2, to=75, by=4), labels=FALSE)
#axis(side=1, at=seq(from=1, to=75, by=4), labels=FALSE)
#axis(side=1, at=seq(from=3, to=75, by=4), labels=FALSE)

#y-axis annotation
mtext(text=expression("Per Cell Respiration Rate (fmol O"[2]*" hr"^-1*")"), 
     side=2, line=1.75, xpd=NA, srt=90)
axis(side=2, 
     at=c(1e-3,1e-2,1e-1,1,1e1,1e2), 
     labels=c(expression("1x10"^-3), expression("1x10"^-2), expression("1x10"^-1), 
              expression("1x10"^0), expression("1x10"^1), expression("1x10"^2)),
     las=1, hadj=0.6, padj=0.40)
axis(side=2, 
     at=c(c(2:9)*1e-3, c(2:9)*1e-2, c(2:9)*1e-1, c(2:9), c(2:9)*1e1), 
     labels=FALSE, tcl=-0.25)

par(ps=8)
text(x=-16., y=2e2, labels="D", cex=1, adj=c(0,0), xpd=NA, font=2)

#third plot
input_file <- paste0(pathDF,"total_O2_consumption_per_ml_per_genus.csv")
input <- read.csv(file=input_file, header=TRUE)
#get unique genus and total abundance
unique_genes <- unique(input$genus)
unique_dates <- unique(input$date)

#reformat to a table with genes and then three days
data_tableA <- data.frame(unique_genes, 
                         as.numeric(input$Total_O2_consumption[input$date == unique_dates[1]]), 
                         as.numeric(input$Total_O2_consumption[input$date == unique_dates[2]]), 
                         as.numeric(input$Total_O2_consumption[input$date == unique_dates[3]]))
colnames(data_tableA) <- c("genus",unique_dates)
data_tableB <- merge(x=data_table3, y=data_tableA, by.x="genus", by.y="genus",
                     sort=FALSE, all=FALSE)
data_tableC <- select(data_tableB, "genus", "Color", "class", "20181030", "20190402", "20190709", "ClassColor")
colnames(data_tableC) <- c("genus", "Color", "class", "X20181030", "X20190402", "X20190709", "ClassColor")

#basic plot
par(fig=c(0.15,0.99,ylim1,ylim1+y_height), new=TRUE)
par(ps=6)
y_data <- merge(x=combined_names, y=data_tableC, by.x="genus", by.y="genus",
                sort=FALSE, all=FALSE)

#combind other and unclassified into one category
other_1 <- 1 - sum(y_data$X20181030)/sum(data_tableC$X20181030)
other_2 <- 1 - sum(y_data$X20190402)/sum(data_tableC$X20190402)
other_3 <- 1 - sum(y_data$X20190709)/sum(data_tableC$X20190709)

unclass_1 <- y_data$X20181030[19]/sum(data_tableC$X20181030)
unclass_2 <- y_data$X20190402[19]/sum(data_tableC$X20190402)
unclass_3 <- y_data$X20190709[19]/sum(data_tableC$X20190709)

other_O2 <- c(other_1+unclass_1, other_2+unclass_2, other_3+unclass_3)

x_seq <- seq(from=2, to=74, by=4)   #1:nrow(y_data)
points2D(x=x_seq[1:18]-1, 
         y=(100*y_data$X20181030/sum(data_tableC$X20181030))[1:18], 
         ylim=c(1e-3,1e2), 
         xlim=c(1,4*length(y_data$X20181030)-2),
         xaxt='n', yaxt='n', pch=15, xlab="", ylab="", 
         col=y_data$ClassColor, log="y", cex=0.4)
points2D(x=x_seq[1:18], 
         y=(100*y_data$X20190402/sum(data_tableC$X20190402))[1:18], 
         xaxt='n', yaxt='n', pch=15, 
         col=y_data$ClassColor, log="y", cex=0.4, add=TRUE)
points2D(x=x_seq[1:18]+1, 
         y=(100*y_data$X20190709/sum(data_tableC$X20190709))[1:18], 
         xaxt='n', yaxt='n', pch=15, 
         col=y_data$ClassColor, log="y", cex=0.4, add=TRUE)
points2D(x=c(4*nrow(y_data)-c(1,2,3)), 
         y=100*other_O2, 
         xaxt='n', yaxt='n', pch=15, 
         col=other_color, log="y", cex=0.4, add=TRUE)
#add grey vertical lines between genus
for (iline in seq(from=0, to=79, by=4)){
  lines(x=c(iline,iline), y=c(0.001,100), type="l", col="grey50", lwd=0.5)
}
#x-axis annotation
#x_labels <- y_data$genus
x_labels <- c("Pelagibacter", 
              "SCGC-AAA076-P13", 
              "MS024-2A", 
              "D2472", 
              "Thioglobus", 
              "Hel1-33-131",    
              "BACL14", 
              "Planktomarina",   #starrred
              "Amylibacter",     #starrred
              "MAG-121220-bin8", 
              "Polaribacter",    #starrred
              "SCGC-AAA160-P02", #starrred
              "HC6-5",           #starrred
              "GCA-002733185",   #starrred
              "UBA7428",         #starrred
              "ASP10-02a",       #starrred
              "UBA4465",         #starrred
              "LFER01",          #starrred
              "Unclassified")
x_labels[19] <- "Others"
x_colors <- c(rep("grey30",7),rep("black",2),"grey30",rep("black",8),"grey30")   #c(y_data$name_color,other_color)
x_ticks <- c(x_seq,78)  #1:(1+nrow(y_data))
text(x=x_ticks[1:19], 
     y=0.3e-3, 
     adj=c(1,1),
     labels=x_labels[1:19], 
     srt=45, 
     cex=1, 
     xpd=NA, 
     col=x_colors)
axis(side=1, at=x_ticks, labels=FALSE)
#y-axis annotation
mtext(text=expression("Total O"[2]*" Consumption (daily %)"),
      side=2, line=1.75, xpd=NA, srt=90)
axis(side=2, 
     at=c(1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2), 
     labels=c(expression("1x10"^-3), expression("1x10"^-2), expression("1x10"^-1),
              expression("1x10"^0), expression("1x10"^1), expression("1x10"^2)),
     las=1, hadj=0.6, padj=0.40)
axis(side=2, 
     at=c(c(2:9)*1e-3, c(2:9)*1e-2, c(2:9)*1e-1, c(2:9)*1e0, c(2:9)*1e1), 
     labels=FALSE, tcl=-0.25)

par(ps=8)
text(x=-16., y=2e2, labels="E", cex=1, adj=c(0,0), xpd=NA, font=2)

#close file and reset plotting parameters to the default
par(op)
dev.off()