
###### Setup (load packages and read in data) ########################

#load all relevant packages
library(adegenet)
library(ggplot2)
library(cowplot)
library(factoextra)
library(readr)
library(tidyr)
library(stringr)
library(data.table)
library(dplyr)

#read in color data, in format produced by Colormesh
guppy_color_measurements <- read.csv("guppy_color_measurements.csv")

#read in x and y coordinates for the locations on the images at which color was sampled, as produced by tpsDig
Color_coords <- read.csv("Sampling_locations.csv")


######### Perform the Discriminant Analysis of Principal Components (DAPC) ###############

#n.da is the number of discriminant functions to use, which defaults to 1 less than number of groups
dapc_output = dapc(guppy_color_measurements[,c(3:9906)], guppy_color_measurements$Line, var.contrib = TRUE, scale = TRUE, n.pca = 10, n.da = 10, var.loadings = TRUE) 

#join coordinates from the DAPC with the group labels
dapc_output$ind.coord_lin_names = cbind(guppy_color_measurements$Line, dapc_output$ind.coord)


######### Start making heatmaps describing the meaning of each discriminant function ##########

#define function for plotting dapc results
fish_dapc_loading_heatmap = function(dapc_object, axis, filename, ylim, X_Y_coords, color.channel = "sum",plot.legend=T){
  #extract variable loadings and names from the dapc results, then convert them into data frames
  variable.loadings <- as.data.frame(dapc_object$var.load[,axis])
  variable.names <- as.data.frame(names(dapc_object$var.load[,axis]))
  #join variable loadings and names
  loadings.and.names = cbind(variable.loadings, variable.names)
  colnames(loadings.and.names) = c("corrs", "point")
  #split loadings into the 4 color channels
  red.loadings.and.names = loadings.and.names[1:2476, ]
  green.loadings.and.names = loadings.and.names[2477:4952, ]
  blue.loadings.and.names = loadings.and.names[4953:7428, ]
  UV.loadings.and.names = loadings.and.names[7429:9904, ]
  #combine just the loadings for all channels into one data frame
  all.color.loadings = as.data.frame(cbind(red.loadings.and.names$corrs, green.loadings.and.names$corrs, blue.loadings.and.names$corrs, UV.loadings.and.names$corrs))
  colnames(all.color.loadings) = c("red", "gre", "blu", "UV")
  #combine coordinates and loadings
  coordinates.and.loadings = cbind(X_Y_coords, all.color.loadings)
  #add column with means for each color channel
  coordinates.and.loadings$mean = rowMeans(coordinates.and.loadings[,c("red","gre","blu", "UV")])
  #determine the color scale for the heatmap for each color channel (and for mean color)
  if(color.channel == "mean"){
    all.color.loadings = cbind(X_Y_coords, coordinates.and.loadings[7])
    low.color = "white"
    high.color = "black"
  } else if(color.channel == "red"){
    all.color.loadings = cbind(X_Y_coords, coordinates.and.loadings$red)
    low.color = "black"
    high.color = "red"
  }else if(color.channel == "green"){
    all.color.loadings = cbind(X_Y_coords, coordinates.and.loadings$gre)
    low.color = "black"
    high.color = "green"
  }else if(color.channel == "blue"){
    all.color.loadings = cbind(X_Y_coords, coordinates.and.loadings$blu)
    low.color = "black"
    high.color = "blue"
  }else {   #default is UV
    all.color.loadings = cbind(X_Y_coords, coordinates.and.loadings$UV)
    low.color = "black"
    high.color = "white"
  }
  names(all.color.loadings) = c("V1","V2", "response")
  
  if(plot.legend){
    #create plot with a legend
    ggplot(all.color.loadings, aes(V1,-V2, color = response)) +
      ggtitle("") +
      geom_point(size=1.5) +
      coord_fixed() +
      scale_color_gradient2(low = low.color, mid = "grey", high = high.color,breaks = c(-0.006,0,0.006),labels = c(-0.006,0,0.006), limits = ylim, name = toupper(color.channel)) +
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            plot.title = element_text(size = 24, hjust = 0.5, vjust=5),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),plot.background=element_blank(),
            legend.position = "right",
            legend.text = element_text(size=10),
            legend.title = element_text(size=12))
  } else {
    #create plot without legend 
    ggplot(all.color.loadings, aes(V1,-V2, color = response)) +
      ggtitle("") +
      geom_point(size=1.5) +
      coord_fixed() +
      scale_color_gradient2(low = low.color, mid = "grey", high = high.color, limits = ylim, name = toupper(color.channel)) +
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            plot.title = element_text(size = 24, hjust = 0.5, vjust=5),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),plot.background=element_blank(),
            legend.position = "none") 
  }
  
}

#create vector of color channel names
cols <- c("red","green","blue","UV")

#use heatmap function to create heatmap of fish, once for each color channel, depicting effect of each discriminant function on color pattern
#heatmap for discriminant function 1
dapc1_fish <- lapply(cols, function(colour){
  fish_dapc_loading_heatmap(dapc_object = dapc_output, axis = 1, filename =  paste(output_folder, "DAPC_ax1_red_corrs.png", sep = ""), ylim = c(-0.006,0.006), X_Y_coords = Color_coords, color.channel = colour, plot.legend = F)
})
dapc1_fish[[1]] <- dapc1_fish[[1]] + ggtitle("DF1")
dapc1_fish_fig <- cowplot::plot_grid(plotlist = dapc1_fish, ncol=1)

#heatmap for discriminant function 2
#include legend
dapc2_fish <- lapply(cols,function(colour){
  fish_dapc_loading_heatmap(dapc_object = dapc_output, axis = 2, filename =  paste(output_folder, "DAPC_ax1_red_corrs.png", sep = ""), ylim = c(-0.006,0.006), X_Y_coords = Color_coords, color.channel = colour, plot.legend = T)
})
dapc2_fish[[1]] <- dapc2_fish[[1]] + ggtitle("DF2")
dapc2_fish_fig <- cowplot::plot_grid(plotlist = dapc2_fish, ncol=1)


###### start making DAPC figure ########

#create dataframe with coordinates for each individual from the DAPC
plot_dd <- data.frame(dapc_output$ind.coord)
#add group (line) labels
plot_dd$IF_line <- dapc_output$grp

#create key mapping colors to line names for the DAPC plot
cols_dapc <- data.frame(line = paste0("Iso-Y", c(10, 6, 8, 9)),
                   col = c("#3CBB75", "#B22222", "#2D708E", "#481567"))

#create vector of group labels to display on DAPC plot
plot_labs <- data.frame(dapc_output$grp.coord)
plot_labs$IF_line <- rownames(plot_labs)
plot_labs$IF_line <- c("Iso-Y10", "Iso-Y6", "Iso-Y8", "Iso-Y9")

#get the DAPC eigenvals
eigenvals <- round((dapc_output$eig/sum(dapc_output$eig)) * 100, 2)

dapc_fig <- ggplot(plot_dd, aes(LD1, LD2, colour = IF_line)) +
  theme_classic() +
  geom_vline(xintercept = 0, alpha = 0.3) +
  geom_hline(yintercept = 0, alpha = 0.3) +
  geom_point(show.legend = F, size = 3, alpha = 0.7) +
  scale_colour_manual(breaks = cols_dapc$line,
                      values = cols_dapc$col) +
  geom_segment(data = plot_dd[plot_dd$IF_line == "Iso-Y6",], aes(x = LD1, xend = plot_labs[plot_labs$IF_line == "Iso-Y6", "LD1"],
                                                              y = LD2, yend = plot_labs[plot_labs$IF_line == "Iso-Y6", "LD2"])) +
  geom_segment(data = plot_dd[plot_dd$IF_line == "Iso-Y8",], aes(x = LD1, xend = plot_labs[plot_labs$IF_line == "Iso-Y8","LD1"],
                                                              y = LD2, yend = plot_labs[plot_labs$IF_line == "Iso-Y8", "LD2"])) +
  geom_segment(data = plot_dd[plot_dd$IF_line == "Iso-Y9",], aes(x = LD1, xend = plot_labs[plot_labs$IF_line == "Iso-Y9", "LD1"],
                                                              y = LD2, yend = plot_labs[plot_labs$IF_line == "Iso-Y9", "LD2"])) +
  geom_segment(data = plot_dd[plot_dd$IF_line == "Iso-Y10",], aes(x = LD1, xend = plot_labs[plot_labs$IF_line == "Iso-Y10", "LD1"],
                                                               y = LD2, yend = plot_labs[plot_labs$IF_line == "Iso-Y10", "LD2"])) +
  geom_label(data = plot_labs, aes(x = LD1, y = LD2, label = IF_line), show.legend = F, size = 8, alpha = 0.8) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.position = "none") +
  scale_x_continuous(breaks = c(-6, 0, 6)) +
  scale_y_continuous(breaks = c(-6, 0, 6)) +
  labs(x = paste0("Discriminant Function 1 (", eigenvals[1], "%)"),
       y = paste0("Discriminant Function 2 (", eigenvals[2], "%)")) +
  stat_ellipse()


############# Create and save combined DAPC and heatmaps figure ##############

empty <- ggplot() + theme_void()

png("dapc_and_heatmaps.png", width = 40, height = 20, units = "cm", res = 400)
plot_grid(dapc_fig, dapc1_fish_fig, dapc2_fish_fig, ncol = 3, rel_widths = c(2.5, 0.75, 1.05))
dev.off()


