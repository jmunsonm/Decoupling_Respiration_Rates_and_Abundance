rm(list = ls())

library(plot3D)
library(dplyr)
library(vioplot)

path0 <- "~/Desktop/Git_upload/"
pathDF <- "~/Desktop/Git_upload/Data Files/"
pathFIG <- "~/Desktop/Git_upload/Figures/"

log_model_Fig1B <- function(Winkler){
  model2 <- lm(log10(Winkler.method) ~ log10(RSG.method), data = Winkler)
  summary2 <- summary(model2)
  slope <- summary2$coefficients[2,1]
  intercept <- summary2$coefficients[1,1]
  model_Rsq <- summary2$r.squared
  
  #add confidence intervals on prediction
  x2_model <- data.frame(10^seq(from=-3, to=1, by=0.1))
  colnames(x2_model) <- "RSG.method"
  y2_model <- 10^(predict(model2, newdata=x2_model, interval="prediction"))
  df <- data.frame(x2_model,y2_model)
  
  #CI on prediction
  x_CI <- c(x2_model[,1],x2_model[nrow(x2_model):1,1])
  y_CI <- c(df[,3],df[nrow(df):1,4])
  df_CI <- data.frame(x_CI, y_CI)
  polygon2D(x=df_CI$x_CI, y=df_CI$y_CI, log="xy",
            xaxt="n", yaxt='n', 
            ylab="", 
            xlab=expression(paste("RSG-based Respiration Rate ("*mu*"mol O"[2]~"L"^-1~"h"^-1*")")),
            xlim=c(0.001, 10), ylim=c(0.001, 10),
            col="grey90", facets=TRUE, add=FALSE)
  
  #add confidence intervals on mean
  x2_model <- data.frame(10^seq(from=-3, to=1, by=0.1))
  colnames(x2_model) <- "RSG.method"
  y2_model <- 10^(predict(model2, newdata=x2_model, interval="confidence"))
  df <- data.frame(x2_model,y2_model)
  
  #CI on mean
  x_CI <- c(x2_model[,1],x2_model[nrow(x2_model):1,1])
  y_CI <- c(df[,3],df[nrow(df):1,4])
  df_CI <- data.frame(x_CI, y_CI)
  polygon2D(x=df_CI$x_CI, y=df_CI$y_CI, log="xy",
            xlim=c(0.001, 10), ylim=c(0.001, 10),
            col="grey", facets=TRUE, add=TRUE)
  
  points2D(x=Winkler$RSG.method, 
           y=Winkler$Winkler.method, 
           add=TRUE,
           col=Winkler$Color, pch=Winkler$Shape, las=1)
  
  lines2D(x=df[,1], y=df[,2], col="black", lwd=1.5, log="y", add=TRUE)
  text(x=0.1, y=0.5e-2, 
       labels=bquote("log"[10]*"(y) =" ~ .(round(intercept,3)) ~ 
                       "+" ~ .(round(slope,3))~"log"[10]*"(x)"),
       xpd=NA, adj=c(0,0.25), srt=0, cex=1.5)
  text(x=0.1, y=0.5e-2, 
       labels=bquote("R"^2 ~ "=" ~ .(round(model_Rsq,3))),
       xpd=NA, adj=c(0,1.25), srt=0, cex=1.5)
  legend(x="topleft", pch=15, 
         legend=c("GoM", "Atlantic Ocean", "Mallorca"),
         col=c("#003f5c", "#7a5195", "#ef5675"),
         y.intersp = 0.75, title="Source", bty="n")
  legend(x="top", pch=c(15,16), legend=c("Open Ocean", "Coastal"), 
         y.intersp = 0.75, title="Ocean Location", bty="n")
  
  #x-axis annotation
  axis(side=1, at=10^c(-3:1), las=1, padj=0.,
       labels=c(0.001, 0.01, 0.1, 1.0, 10.))
  axis(side=1, at=(c((2:9)*1e-3,(2:9)*1e-2, (2:9)*1e-1, (2:9)*1e0, (2:9)*1e1)),
       tcl=-0.25, labels=FALSE)
  
  #y-axis annotation
  mtext(text=expression(paste("Winkler-based Respiration Rate")), 
        xpd=NA, side=2, line=2.8, cex=1.1)
  mtext(text=expression(paste("("*mu*"mol O"[2]~"L"^-1~"h"^-1*")")),
        xpd=NA, side=2, line=1.75, cex=1.1)
  axis(side=2, at=10^c(-3:1), las=1, hadj=1.05,
       labels=c(0.001, 0.01, 0.1, 1.0, 10.))
  axis(side=2, at=(c((2:9)*1e-3,(2:9)*1e-2, (2:9)*1e-1, (2:9)*1e0, (2:9)*1e1)),
       tcl=-0.25, labels=FALSE)
  
  par(font=2, ps=12)
  y_tab <- 10^(1+0.075*(1+3))
  text(x=1e-3, y=y_tab, labels="B", xpd=NA, pos=2, offset=4)
  par(font=1, ps=7)
}

linear_model_Fig1B <- function(Winkler, thru_origin=FALSE){
  x <- Winkler$RSG.method
  y <- Winkler$Winkler.method
  
  if(thru_origin){
    model1 <- lm(Winkler.method ~ 0+RSG.method, data = Winkler)
    summary1 <- summary(model1)
    slope <- summary1$coefficients[1,1]
    intercept <- 0
    model_Rsq <- summary1$r.squared
  } else {
  model1 <- lm(Winkler.method ~ RSG.method, data = Winkler)
  summary1 <- summary(model1)
  slope <- summary1$coefficients[2,1]
  intercept <- summary1$coefficients[1,1]
  model_Rsq <- summary1$r.squared
  }
  #add confidence intervals on prediction
  x1_model <- data.frame((0:300)/100)
  colnames(x1_model) <- "RSG.method"
  y1_model <- predict(model1, newdata=x1_model, interval="prediction")
  df <- data.frame(x1_model,y1_model)
  
  #CI on prediction
  x_CI <- c(x1_model[,1],x1_model[nrow(x1_model):1,1])
  y_CI <- c(df[,3],df[nrow(df):1,4])
  df_CI <- data.frame(x_CI, y_CI)
  polygon2D(x=df_CI$x_CI, y=df_CI$y_CI, 
            xlim=c(0.00, 3), ylim=c(0.00, 3), xaxt='n', yaxt='n',
            xlab=expression(paste("RSG-based Respiration Rate ("*mu*"mol O"[2]~"L"^-1~"h"^-1*")")),
            ylab="",
            col="grey90", facets=TRUE, add=FALSE)
  axis(side=1, at=0:3, labels=TRUE)
  axis(side=1, at=c(0.2,0.4,0.6,0.8, 1+c(0.2,0.4,0.6,0.8), 2+c(0.2,0.4,0.6,0.8)),
       labels=FALSE, tcl=-0.25)
  axis(side=2, at=0:3, labels=TRUE, las=1, hadj=1.5)
  axis(side=2, at=c(0.2,0.4,0.6,0.8, 1+c(0.2,0.4,0.6,0.8), 2+c(0.2,0.4,0.6,0.8)),
       labels=FALSE, tcl=-0.25)
  
  #add confidence intervals on mean
  y2_model <- predict(model1, newdata=x1_model, interval="confidence")
  df <- data.frame(x1_model,y2_model)
  
  #CI on mean
  x_CI <- c(x1_model[,1],x1_model[nrow(x1_model):1,1])
  y_CI <- c(df[,3],df[nrow(df):1,4])
  df_CI <- data.frame(x_CI, y_CI)
  polygon2D(x=df_CI$x_CI, y=df_CI$y_CI, log="xy",
            xlim=c(0.001, 10), ylim=c(0.001, 10),
            col="grey", facets=TRUE, add=TRUE)
  
  points2D(x=x, 
           y=y, 
           xlim=c(0.00, 3), ylim=c(0.00, 3),
           add=TRUE,
           xlab=expression(paste("RSG-based Respiration Rate ("*mu*"mol O"[2]~"L"^-1~"h"^-1*")")),
           ylab="",
           col=Winkler$Color, pch=Winkler$Shape, las=1)
  mtext(text="Winkler-based Respiration Rate", side=2, line=3, cex=1.2)
  mtext(text=expression("("*mu*"mol O"[2]~"L"^-1~"h"^-1*")"), side=2, line=1.5, cex=1.2)
  lines2D(x=df[,1], y=df[,2], col="red", lwd=2, add=TRUE)
  if(thru_origin){
    text(x=0.25, y=1.75, 
         labels=bquote("y =" ~ .(round(slope,3))~"*x"),
         xpd=NA, adj=c(0,0), srt=0, cex=1.5)
  } else {
    text(x=0.25, y=1.75, 
         #labels=bquote("y =" ~ .(round(slope,3))~"*x"),
         labels=bquote("y =" ~ .(round(intercept,3))~"+"~ .(round(slope,3))~"*x"),
         xpd=NA, adj=c(0,0), srt=0, cex=1.5)
  }
  text(x=0.25, y=1.75, 
       labels=bquote("R"^2 ~ "=" ~ .(round(model_Rsq,3))),
       xpd=NA, adj=c(0,1), srt=0, cex=1.5)
  lines2D(x=(0:30)/10, y=(0:30)/10, col="black", lty=2, lwd=2, add=TRUE)
  legend(x="topleft", pch=15, 
         legend=c("GoM", "Atlantic Ocean", "Mallorca"), 
         col=c("#003f5c", "#7a5195", "#ef5675"),
         y.intersp = 0.75, title="Source", bty="n")
  legend(x="top", pch=c(15,16), legend=c("Open Ocean", "Coastal"), 
         y.intersp = 0.75, title="Ocean Location", bty="n")
  
  par(font=2, ps=12)
  text(x=0, y=3.2, labels="B", xpd=NA, pos=2, offset=4)
  par(font=1, ps=7)
  
  #add inset figure
  par(fig=c(0.779,0.974,0.815,0.88), mar=c(0,0,0,0), bg="white",new=TRUE)
  points2D(x=Winkler$RSG.method[1:8],
           y=Winkler$Winkler.method[1:8], 
           xlim=c(0,0.12), ylim=c(0,0.10),
           xlab="", ylab="", cex.axis=1.0, 
           col=Winkler$Color[1:8],
           pch=Winkler$Shape[1:8],
           add=FALSE)
  lines2D(x=df[,1], y=df[,2], col="red", lwd=2, add=TRUE)
  lines2D(x=c(0,0.12), y=c(0,0.12), col="black", lty=2, add=TRUE)
}

#open output plot file
plot_out <- paste0(pathFIG,"Figure 1.pdf")
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

par(mfg=c(1,1))

#data input and cleaning
input_file <- paste0(pathDF,"rsg_reg_data_clean.csv")
rsq_reg <- read.csv(file=input_file, header=TRUE)
keep <- !is.na(rsq_reg[,1])
rsq_reg <- rsq_reg[keep,]

#set plotting shape
shape <- rep(15, nrow(rsq_reg))
change <- rsq_reg$Temperature == 22
shape[change] <- 16

#set plotting colors
base_colors <- c("#F16745", "#FFC65D", "#000099", "#7BC8A4", "#4CC3D9", "#93648D", "#FF33CC")
#"Shewanella"  "Pseudomonas" "Ruegaria"    "Roseobacter" "Oligotropha" "Rhodoluna"   "Vibrio"
color <- rep(base_colors[1], nrow(rsq_reg))
color[rsq_reg$Species == "Pseudomonas"] <- base_colors[2]
color[rsq_reg$Species == "Ruegaria"] <- base_colors[3]
color[rsq_reg$Species == "Roseobacter"] <- base_colors[4]
color[rsq_reg$Species == "Oligotropha"] <- base_colors[5]
color[rsq_reg$Species == "Rhodoluna"] <- base_colors[6]
color[rsq_reg$Species == "Vibrio"] <- base_colors[7]

x <- rsq_reg$rsg_gm
y <- rsq_reg$resp_complexcorrection  ############  NAME CHANGE???
xy_df <- data.frame(x,y)

#linear regression model
rg_model <- lm(log10(y)~x, data=xy_df)
rg_summary <- summary(rg_model)
slope <- rg_summary$coefficients[2,1]
intercept <- rg_summary$coefficients[1,1]
model_Rsq <- rg_summary$r.squared

x1_model <- data.frame(x=seq(from=0, to=35, by=1))

#predicted value and confidence interval on the prediction
y2_model <- 10^(predict(rg_model, newdata=x1_model, interval="prediction"))
df <- data.frame(x1_model,y2_model)
x_CI <- c(x1_model[,1],x1_model[nrow(x1_model):1,1])
y_CI <- c(df[,3],df[nrow(df):1,4])
df_CI <- data.frame(x_CI, y_CI)
polygon2D(x=df_CI$x_CI, y=df_CI$y_CI, log="y",
          col="grey90", facets=TRUE,
          xlim=c(0,35), ylim=c(1e-3,1e1),
          xlab="Normalized RSG Fluorescence Intensity",
          ylab="",
          xaxt="n", yaxt='n')

#predicted value and confidence interval on the mean
y2_model <- 10^(predict(rg_model, newdata=x1_model, interval="confidence"))
df <- data.frame(x1_model,y2_model)

#CI on mean
x_CI <- c(x1_model[,1],x1_model[nrow(x1_model):1,1])
y_CI <- c(df[,3],df[nrow(df):1,4])
df_CI <- data.frame(x_CI, y_CI)
polygon2D(x=df_CI$x_CI, y=df_CI$y_CI, log="y",
          col="grey", facets=TRUE, add=TRUE)

lines2D(x=df[,1], y=df[,2], col="black", lwd=1.5, log="y", add=TRUE)

#data points themselves
points2D(x=x, y=y, log="y", 
         col=color, pch=shape, add=TRUE)

#par(font=3)
legend(x="bottomright", bty="n", title="Species",
       legend=c("Shewanella loihica", 
                "Pseudomonas stutzeri", 
                "Ruegaria pomeroyi",
                "Roseobacter denitrificans", 
                "Oligotropha carboxidovorans", 
                "Rhodoluna lacicola", 
                "Vibrio alginolyticus"),
       pch=15, text.font=3,
       col=base_colors,  #c("#F16745", "#FFC65D", "#7BC8A4", "#4CC3D9", "#93648D", "#000099"),     #    "#404040"), 
       y.intersp = 0.65)
legend(x="bottom", bty="n",
       legend=c("14C", "22C"), pch=c(15,16), col="black", title="Temperature", 
       y.intersp = 0.65)
par(font=1)

#x-axis annotation
x_ticks <- seq(from=0, to=35, by=5)
axis(side=1, at=x_ticks, labels=TRUE)
axis(side=1, at=c(1:4,6:9,11:14,16:19,21:24,26:29,31:34), tcl=-0.25, labels=FALSE)

#y-axis annotation
par(mgp=c(3,0.6,0))
mtext(text=expression(paste("Respiration Rate (fmol h"^-1~"cell"^-1*")")), 
      xpd=NA, side=2, line=2.6, cex=1.1)
axis(side=2, at=10^c(-3:1), las=1,
   labels=c(0.001, 0.01, 0.1, 1.0, 10))
axis(side=2, at=c(c(2:9)*1e-3, c(2:9)*1e-2, c(2:9)*1e-1, c(2:9)*1e0, c(2:9)*1e1), 
     labels=FALSE, tcl=-0.25)

par(mgp=c(1.75,0.5,0))

#regression line information
y_reg <- 10^(-3+0.95*(1+3))
text(x=1, y=y_reg, 
     labels=bquote("log"[10]*"(y) =" ~ .(round(intercept,3)) ~ 
                      "+" ~ .(round(slope,3))*"(x)"),
     xpd=NA, adj=c(0,0.25), srt=0, cex=1.5)
text(x=1, y=y_reg, 
     labels=bquote("R"^2 ~ "=" ~ .(round(model_Rsq,3))),
     xpd=NA, adj=c(0,1.25), srt=0, cex=1.5)

par(font=2, ps=12)
y_tab <- 10^(1+0.075*(1+3))
text(x=0, y=y_tab, labels="A", xpd=NA, pos=2, offset=4)
par(font=1, ps=7)


#####################################################################
#RSG O2 vs Winker plot
input_file <- paste0(pathDF,"Winkler_RSG_validation_measurements.csv")
Winkler <- read.csv(file=input_file, header=TRUE)
par(mfg=c(1,2))
change <- (Winkler$Color == "#ffa600")
Winkler$Color[change] <- "#ef5675"
#linear_model_Fig1B(Winkler, thru_origin=TRUE)
log_model_Fig1B(Winkler)
par(op)
dev.off()
