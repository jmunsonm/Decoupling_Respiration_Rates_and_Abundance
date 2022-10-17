
plot_logxy <- function(plot_data,
                       x_lim, x_main_label, major_x_tic_location, major_x_tic_labels, minor_x_tic_location,
                       y_lim, y_main_label, major_y_tic_location, major_y_tic_labels, minor_y_tic_location,
                       reg_eqn_xy, pt_label, label_points_xy, tab_y, tab_label){
  points2D(x=plot_data$x, 
           y=plot_data$y, 
           xlab=x_main_label, 
           ylab="",
           col=plot_data$color, 
           log="xy", 
           xlim=x_lim, xaxt='n', 
           ylim=y_lim, yaxt='n',
           pch=plot_data$shape, xpd=NA)
  
  #x-axis annotation
  axis(side=1, at=major_x_tic_location, labels=major_x_tic_labels)
  axis(side=1, at=minor_x_tic_location, labels=FALSE, tcl=-0.25)
  
  #y-axis annotation
  mtext(text=y_main_label, side=2, line=2.6, xpd=NA, cex=1.1)
  axis(side=2, at=major_y_tic_location, las=1, labels=major_y_tic_labels)
  axis(side=2, at=minor_y_tic_location, labels=FALSE, tcl=-0.25)
  
  #add average line
  median_y <- median(plot_data$y)
  x <- x_lim
  y <- rep(median_y,2)
  lines2D(x=x, y=y, col="black", lwd=1.5, lty=2, add=TRUE)
  
  #add regression line
  if(!is.null(reg_eqn_xy)){
    log_x <- log10(plot_data$x)
    log_y <- log10(plot_data$y)
    plot_data$log_x <- log_x
    plot_data$log_y <- log_y
    my_model <- lm(log_y ~ log_x, data = plot_data)
    intercept <- my_model$coefficients[1]
    slope <- my_model$coefficients[2]
    model_Rsq <- summary(my_model)$r.squared
    x <- x_lim
    y <- 10^(intercept + log10(x)*slope)
    lines2D(x=x, y=y, col="black", lwd=1.5, add=TRUE)
    
    #regression line information
    if(slope >= 0){
      text(x=reg_eqn_xy[1], y=reg_eqn_xy[2], 
           labels=bquote("log"[10]*"(y) =" ~ .(round(intercept,3)) ~ 
                           "+" ~ .(round(abs(slope),3))*"*log"[10]*"(x)"),
           xpd=NA, adj=c(0,0.25), srt=0, cex=1.5)
    } else{
      text(x=reg_eqn_xy[1], y=reg_eqn_xy[2], 
           labels=bquote("log"[10]*"(y) =" ~ .(round(intercept,3)) ~ 
                           "-" ~ .(round(abs(slope),3))*"*log"[10]*"(x)"),
           xpd=NA, adj=c(0,0.25), srt=0, cex=1.5)
    }
    text(x=reg_eqn_xy[1], y=reg_eqn_xy[2], 
         labels=bquote("R"^2 ~ "=" ~ .(round(model_Rsq,3))),
         xpd=NA, adj=c(0,1.25), srt=0, cex=1.5)
  }
  
  if(pt_label){
    label_points(data_x=plot_data$x, 
                 data_y=plot_data$y, 
                 data_genus=plot_data$genus, 
                 label_x=label_points_xy[1,],
                 label_y=label_points_xy[2,])
  }
  par(font=2, ps=12)
  text(x=x_lim[1], y=tab_y, labels=tab_label, xpd=NA, pos=2, offset=4)
  par(font=1, ps=7)
}

plot_logy <- function(plot_data,
                      x_lim, x_main_label, major_x_tic_location, major_x_tic_labels, minor_x_tic_location,
                      y_lim, y_main_label, major_y_tic_location, major_y_tic_labels, minor_y_tic_location,
                      reg_eqn_xy, pt_label, label_points_xy, tab_y, tab_label){
  points2D(x=plot_data$x,
           y=plot_data$y, 
           xaxt="n", yaxt='n', 
           xlab=x_main_label, 
           ylab="", 
           xlim=x_lim,
           col=plot_data$color, 
           log="y", ylim=y_lim, 
           pch=plot_data$shape, xpd=NA)
  
  #add regression line
  plot_data$log_y <- log10(plot_data$y)
  X16S_model <- lm(log_y ~ x, data = plot_data)
  intercept <- X16S_model$coefficients[1]
  slope <- X16S_model$coefficients[2]
  model_Rsq <- summary(X16S_model)$r.squared
  x <- x_lim
  y <- 10^(intercept + x*slope)
  lines2D(x=x, y=y, col="black", lwd=1.5, add=TRUE)
  
  #add average line
  median_X16S <- median(plot_data$y)
  x <- x_lim
  y <- rep(median_X16S,2)
  lines2D(x=x, y=y, col="black", lwd=1.5, lty=2, add=TRUE)
  
  #x-axis annotation
  x_ticks <- major_x_tic_location
  axis(side=1, at=major_x_tic_location, labels=major_x_tic_labels)
  axis(side=1, at=minor_x_tic_location, labels=FALSE, tcl=-0.25)
  
  #y-axis annotation
  mtext(text=y_main_label, side=2, xpd=NA, line=2.6, cex=1.1)
  axis(side=2, at=major_y_tic_location, labels=major_y_tic_labels, las=1)
  axis(side=2, at=c(c(2:9)*1e1, c(2:9)*1e2, c(2:9)*1e3, c(2:9)*1e4), 
       labels=FALSE, tcl=-0.25)
  
  #regression line information
  text(x=reg_eqn_xy[1], y=reg_eqn_xy[2], 
       labels=bquote("log"[10]*"(y) =" ~ .(round(intercept,3)) ~ "+" ~ .(round(slope,3))*"*(x)"),
       xpd=NA, adj=c(0,0.25), srt=0, cex=1.5)
  text(x=reg_eqn_xy[1], y=reg_eqn_xy[2], 
       labels=bquote("R"^2 ~ "=" ~ .(round(model_Rsq,3))),
       xpd=NA, adj=c(0,1.25), srt=0, cex=1.5)
  
  if(pt_label){
    label_points(data_x=plot_data$x, 
                 data_y=plot_data$y, 
                 data_genus=plot_data$genus, 
                 label_x=label_points_xy[1,],
                 label_y=label_points_xy[2,])
  }
  par(font=2, ps=12)
  text(x=1.5e-3, y=tab_y, labels=tab_label, xpd=NA, pos=2, offset=4)
  par(font=1, ps=7)
}

plot_logx <- function(plot_data,
                      x_lim, x_main_label, major_x_tic_location, major_x_tic_labels, minor_x_tic_location,
                      y_lim, y_main_label, major_y_tic_location, major_y_tic_labels, minor_y_tic_location,
                      reg_eqn_xy, pt_label, label_points_xy, tab_y, tab_label){
  points2D(x=plot_data$x, 
           y=plot_data$y, 
           xaxt="n", yaxt='n', 
           xlab=x_main_label, 
           ylab="",
           xlim=x_lim, ylim=y_lim,
           col=plot_data$color, 
           log="x", cex.axis=0.7, pch=plot_data$shape, xpd=NA)
  
  #x-axis annotation
  axis(side=1, at=major_x_tic_location, las=1, labels=major_x_tic_labels)
  axis(side=1, at=minor_x_tic_location, labels=FALSE, tcl=-0.25)
  
  #y-axis annotation
  mtext(text=y_main_label, xpd=NA, side=2, line=2.6, cex=1.1)
  axis(side=2, at=major_y_tic_location, las=1, labels=major_y_tic_labels)
  axis(side=2, at=minor_y_tic_location, labels=FALSE, tcl=-0.25)
  
  #add average line
  median_y <- median(plot_data$y)
  x <- x_lim
  y <- rep(median_y,2)
  lines2D(x=x, y=y, col="black", lwd=1.5, lty=2, add=TRUE)
  
  #add regression line
  plot_data$log_x <- log10(plot_data$x)
  my_model <- lm(y ~ log_x, data = plot_data)
  intercept <- my_model$coefficients[1]
  slope <- my_model$coefficients[2]
  model_Rsq <- summary(my_model)$r.squared
  x <- x_lim
  y <- intercept + log10(x)*slope
  lines2D(x=x, y=y, col="black", lwd=1.5, add=TRUE)
  
  #regression line information
  if (slope >= 0){
    text(x=reg_eqn_xy[1], y=reg_eqn_xy[2], 
         labels=bquote("y ="~.(round(intercept,3))~ 
                         "+"~.(round(slope,3))*"*log"[10]*"(x)"),
         xpd=NA, adj=c(0,0.25), srt=0, cex=1.5)
  } else {
    text(x=reg_eqn_xy[1], y=reg_eqn_xy[2], 
         labels=bquote("y ="~.(round(intercept,3))~ 
                         .(round(slope,3))*"*log"[10]*"(x)"),
         xpd=NA, adj=c(0,0.25), srt=0, cex=1.5)
  }
  text(x=reg_eqn_xy[1], y=reg_eqn_xy[2], 
       labels=bquote("R"^2~"="~.(round(model_Rsq,3))),
       xpd=NA, adj=c(0,1.25), srt=0, cex=1.5)
  
  if(pt_label){
    label_points(data_x=plot_data$x, 
                 data_y=plot_data$y, 
                 data_genus=plot_data$genus, 
                 label_x=label_points_xy[1,],
                 label_y=label_points_xy[2,])
  }
  
  par(font=2, ps=12)
  text(x=x_lim[1], y=tab_y, labels=tab_label, xpd=NA, pos=2, offset=4)
  par(font=1, ps=7)
}

find_adj <- function(points_x, points_y, label_x, label_y){
  npts <- length(points_x)
  #number of  points to the left, right, above, and below the label, respectively
  lft <- sum(points_x <= label_x)
  rht <- sum(points_x >= label_x)
  abv <- sum(points_y >= label_y)
  blw <- sum(points_y <= label_y)
  
  if(lft == npts & blw == npts) {
    my_adj <- c(-0.,-0.)                #left lower 
  } else if (lft == npts & blw != npts & abv != npts) {
    my_adj <- c(-0.,0.5)                #left center
  } else if (lft == npts & abv == npts){ 
    my_adj <- c(-0.,1.)                 #left upper 
  } else if(rht == npts & blw == npts){ 
    my_adj <- c(1,-0.0)                 #right lower
  } else if(rht == npts & blw != npts & abv != npts){ 
    my_adj <- c(1.,0.5)                 #right center 
  } else if(rht == npts & abv == npts){ 
    my_adj <- c(1.,1.)                  #right upper 
  } else if(lft != npts & rht != npts & blw == npts){ 
    my_adj <- c(0.5,-0.3)               #center lower
  } else if(lft != npts & rht != npts & abv == npts) {
    my_adj <- c(0.5,1.01)               #center upper
  } else {
    my_adj <- c(0.5,0.5)                #pure center
  }
  return(my_adj)
}

label_points <- function(data_x, data_y, data_genus, label_x, label_y){
  #lines between Pelagibacterales points and its label
  #Pelagibacterales	Pelagibacter
  #Pelagibacterales	Pelagibacter_A
  this_genus <- (data_genus == "Pelagibacter") | (data_genus == "Pelagibacter_A")
  points_x <- data_x[this_genus]
  points_y <- data_y[this_genus]
  if(length(points_x >= 1)){
    text(x=label_x[1], y=label_y[1],
         labels="Pelagibacterales", 
         adj=find_adj(points_x, points_y, label_x[1], label_y[1]),  #c(-0.03,0.5), 
         xpd=NA, cex=1.5)
    for (i in 1:length(points_y)){
      lines2D(x=c(points_x[i],label_x[1]), y=c(points_y[i],label_y[1]), 
              col="black", lwd=0.5, add=TRUE)
    }
  }
  #lines between SAR 86 points and its label
  #SAR 86	SCGC-AAA076-P13
  #SAR 86	D2472
  this_genus <- (data_genus == "SCGC-AAA076-P13") | (data_genus == "D2472")
  points_x <- data_x[this_genus]
  points_y <- data_y[this_genus]
  if(length(points_x >= 1)){
    text(x=label_x[2], y=label_y[2],
         labels="SAR 86",  
         adj=find_adj(points_x, points_y, label_x[2], label_y[2]),  #c(-0.03,0.5),
         xpd=NA, cex=1.5)
    for (i in 1:length(points_y)){
      lines2D(x=c(points_x[i],label_x[2]), y=c(points_y[i],label_y[2]), 
              col="black", lwd=0.5, add=TRUE)
    }
  }
  #lines between Poseidoniales points and its label
  #Poseidoniales	MGIIb-O2
  #Poseidoniales	MGIIa-L1
  #Poseidoniales	Poseidonia
  #Poseidoniales	MGIIa-L2
  this_genus <- (data_genus == "MGIIb-O2") |
    (data_genus == "MGIIa-L1") |
    (data_genus == "Poseidonia") |
    (data_genus == "MGIIa-K1") |
    (data_genus == "MGIIa-L2") 
  points_x <- data_x[this_genus]
  points_y <- data_y[this_genus]
  if(length(points_x) >= 1){
    text(x=label_x[3], y=label_y[3],
         labels="Poseidoniales",  
         adj=find_adj(points_x, points_y, label_x[3], label_y[3]),  #c(-0.03,0.5),
         xpd=NA, cex=1.5)
    for (i in 1:length(points_y)){
      lines2D(x=c(points_x[i],label_x[3]), y=c(points_y[i],label_y[3]), 
              col="black", lwd=0.5, add=TRUE)
    }
  }
  #lines between Rhodobacterales points and its label
  #Rhodobacterales	Amylibacter
  #Rhodobacterales	Planktomarina
  this_genus <- (data_genus == "Planktomarina") | (data_genus == "Amylibacter")
  points_x <- data_x[this_genus]
  points_y <- data_y[this_genus]
  if(length(points_x) >= 1){
    text(x=label_x[4], y=label_y[4],
         labels="Rhodobacterales", 
         adj=find_adj(points_x, points_y, label_x[4], label_y[4]),  #c(-0.03,0.5),
         xpd=NA, cex=1.5)
    for (i in 1:length(points_y)){
      lines2D(x=c(points_x[i],label_x[4]), y=c(points_y[i],label_y[4]), 
              col="black", lwd=0.5, add=TRUE)
    }
  }
}
