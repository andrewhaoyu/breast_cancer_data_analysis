theme_Publication <- function(base_size=18) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, )
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.1), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = 16),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size=15), 
            axis.line = element_line(colour="black",linewidth=2),
            axis.ticks = element_line(),
            # panel.grid.major = element_line(colour="#f0f0f0"),
            # panel.grid.minor = element_line(colour="#f0f0f0"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            #legend.position = "bottom",
            #legend.direction = "horizontal",
            #legend.key.size= unit(0.2, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="bold.italic", size =18),
            #legend.text = element_text(face ="bold"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


scale_fill_Publication <- function(...) {
  library(scales)
  discrete_scale("fill", "Publication", manual_pal(values = c(
    "#386cb0", "#EF7E3D", "#ffd558", "#7fc97f", "#ef3b2c",
    "#662506", "#a6cee3", "#fb9a99", "#984ea3", "#ffff33"
  )), ...)
}

scale_colour_Publication <- function(...) {
  library(scales)
  discrete_scale("colour", "Publication", manual_pal(values = c(
    "#386cb0", "#EF7E3D", "#ffd558", "#7fc97f", "#ef3b2c",
    "#662506", "#a6cee3", "#fb9a99", "#984ea3", "#ffff33"
  )), ...)
}

scale_fill_rainbow <- function(...) {
  library(scales)
  discrete_scale("fill", "Publication", manual_pal(values = c(
    "#5EBD3E", "#FFB900", "#F78200", "#E23838", "#973999", "#009cdf"
  )), ...)
}
