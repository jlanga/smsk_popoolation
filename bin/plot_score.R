#!/usr/bin/env Rscript
require(data.table)
require(dplyr)
require(ggplot2)

colors <- c("darkblue", "orange")

# Function to read the files with the Tajma D's calculations
read_scores <- function(filename){
    
    vs <- filename %>%
    fread(  # Read file and select the columns 
        input = .,
        header = FALSE,
        sep = "\t",
        dec= ".",
        na.strings = "na",
        stringsAsFactors = TRUE,
        select = c(1, 2, 3),
        colClasses = c("character", "numeric", "numeric")
    ) %>%
    select(  # Rename columns
        chromosome = V1,
        position = V2,
        score = V3
    ) %>%
    mutate(  # Convert factor to string 
        chromosome = as.character(chromosome)
    ) %>% 
    arrange(
        chromosome, position
    ) %>% 
    mutate(  # Alternate zeros and ones for the color column
         color = chromosome %>% as.numeric() %% 2
    ) %>%
    mutate(  # Assign the color
         color = factor(
             x = color,
             labels =  c("darkblue", "orange")
         )
    )
    
    vs <- as.data.frame(vs)
    
    return(vs)
}

#Function to do all the plotting
plot_score <- function(data, fileout, hlines=NULL){
    
    data$position <- rownames(data) %>% as.numeric
        
    # Get the labels for the X axis (chromosomes)
    plot_labels <- levels(data$chromosome)
    
    # Compute the tick position
    ticks <- data %>%
        group_by(chromosome) %>%
        summarise(pos = mean(position)) # Compute the placing in the x axis of the tick
    
    
    # Generate the plot: x is the number of elements, y: tajima's D, as color use data$color
    Score_plot <-  ggplot( data=data, aes(x= position, y=score), colour=color) +
        # X axis: title, use breaks and as labels the chromosome names
        scale_x_continuous( "Genomic position" , breaks = ticks$pos , labels = plot_labels ) +
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



