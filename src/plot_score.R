#!/usr/bin/env Rscript
require(argparse)

color_palette <- c("darkblue", "orange")

parser <- ArgumentParser()

parser$add_argument(
    "-i", "--input",
    metavar="FILE",
    type="character",
    help="Input file, with the form CHROMOSOME TAB POSITION TAB SCORE"
)

parser$add_argument(
    "-o", "--output",
    metavar="FILE",
    type="character",
    help="Output file. Format will be inferred from extension."
)

parser$add_argument(
    "-n", "--normalize",
    action="store_true",
    help="Perform normalization of values (value - mean)/std_dev",
    default=FALSE
)

args <- parser$parse_args()


require(tidyverse)
require(ggplot2)



# Function to read the files with the Tajma D's calculations
read_scores <- function(filename){
    
    data <- filename %>%
      read_tsv(  # Read file and select the columns 
          file = .,
          col_names = c("chromosome", "position", "score"),
          col_types = "cnn",
          na = c("na")
      ) %>%
      mutate(  # Convert chromosome to a factor
        chromosome = factor(
          chromosome, 
          levels = chromosome %>% unique()
        )
      ) %>%
      mutate(  # Give alternating numbers
        color = chromosome %>%
            as.numeric() %% 2 + 1
      ) %>%
      mutate( # Give a color
        color = color_palette[color]
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
    plot_labels <- data$chromosome %>% unique()
      
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
        # Y axis: title23 args <- parser$parse_args()
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




#args <- commandArgs(trailingOnly = TRUE)

filein  <- args$input
fileout <- args$output
normalize <- args$normalize

data <- read_scores(filein)

if( normalize ){
    data %>%
        mutate(
            score = (score - mean(score, na.rm= TRUE)) / sd(score, na.rm = TRUE)
        ) %>%
    plot_score(., fileout)
}else{
    plot_score(data, fileout)
}



