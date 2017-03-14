#!/usr/bin/env Rscript
require(tidyverse)
require(ggplot2)

color_palette <- c("darkblue", "orange")

#filename <- "results/tajimad/plot/pop1.tsv"


# Function to read the files with the Tajma D's calculations
read_scores <- function(filename){
    
    data <- filename %>%
    read_tsv(  # Read file and select the columns 
        file = .,
        col_names = c("chromosome", "position", "score"),
        col_types = "cnn",
        na = c("na")
    ) %>%
    mutate(  # Give alternating colors
      color = chromosome %>%
        as.factor() %>%
        as.numeric() %% 2 %>%
        color_palette[.]
    )
    
    data <- as.data.frame(data)
    
    return(data)
}

#Function to do all the plotting
plot_score <- function(data, fileout, hlines=NULL){
    
    data <- data %>% mutate(
      x_position  = row_number()
    )
        
    # Get the labels for the X axis (chromosomes)
    plot_labels <- data %>% 
      group_by(chromosome) %>%
      summarise()
    
    # Compute the tick position
    ticks <- data %>%
        group_by(chromosome) %>%
        summarise(pos = mean(x_position)) # Compute the placing in the x axis of the tick
    
    
    # Generate the plot: x is the number of elements, y: tajima's D, as color use data$color
    Score_plot <-  ggplot(
          data=data, 
          aes(
            x= x_position,
            y= score,
            color= color
          )
        ) +
        # X axis: title, use breaks and as labels the chromosome names
        scale_x_continuous(
          name = "Genomic position",
          breaks = ticks$pos,
          labels = plot_labels ) +
        # Y axis: title
        scale_y_continuous( "Score" ) +
        # Plot points with size 1.5, alpha 0.5 and coloured
        geom_point( aes( alpha = 0.5 , colour = color ) , size = 1.5, na.rm = TRUE ) +
        # Use colors darkblue and orange
        scale_color_manual( values = c( "darkblue", "orange" ) ) +
        # Use gray background
        theme_gray() +
        # Main title
        theme( legend.position = "none" )
    # If supplied constants (horizontal lines) plot them
    if( length( hlines ) > 0 ){
        Score_plot <- Dplot + geom_hline( yintercept = hlines   ,
                                          colour     = "black"  ,
                                          linetype   = "dashed" )
    }
    Score_plot
    ggsave(fileout, width = 293, height = 219, units = "mm")

}




args <- commandArgs(trailingOnly = TRUE)

action  <- args[1]
filein  <- args[2]
fileout <- args[3]

if(length(args) != 3){
    stop("Incorrect number of files. Maybe missing action")
}

data <- read_scores(filein)

if( action == "z"){
    data %>%
        mutate(
            score = (score - mean(score, na.rm= TRUE)) / sd(score, na.rm = TRUE)
        ) %>%
    plot_score( . , fileout )
}else if(action == "none"){
    plot_score( data , fileout )
}else{
    stop("Incorrect action")
}



