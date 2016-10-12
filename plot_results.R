


library(readr)
library(dplyr)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)

model_map <- c("noselection" = "Neutral",
    "empiricalvalues" = "EmpiricalValue",
    "purifyingchem" = "ChemicalAffinity",
    "homoresidue" = "Empirical\nHomoresidue")


plot.fitness.factors <- function(df) {
    
    df <- df %>% 
            group_by(model) %>% 
            mutate(flevel=as.numeric(factor(fitness))) %>% 
            ungroup() %>% 
            mutate(flevel=factor(flevel))
    ggplot(df, aes(model, omega, color=flevel, shape=generation)) +
      theme_bw() +
      geom_point(position = position_jitter(w = 0.25), size=3) +
      scale_colour_brewer(type='qual', palette="Paired", name="Fitness Strength", breaks=c(3, 2, 1),
                          labels=c("high", "middle", "low")) +
      scale_x_discrete( labels=model_map)
      
}

plot.fitness.factors.log <- function(df) {
    # same as plot.fitness.factors(), but used log10 y-axis.
    # doesn't just use scale_y_log10() becaue that results in
    # uninformative axis labels.  Instead
    # breaks = round(seq(min(df$omega), max(df$omega), by=10),1)
    breaks = c(0, 0.5, 1, 2, 10, 20, 30, 40, 50)
    plot.fitness.factors(df) + coord_trans(y="log10") + 
    scale_y_continuous(breaks = breaks)
      
}


df <- read_csv('build/results.csv')
df <- df %>% filter(complete.cases(.)) %>%
       #mutate(omega=replace(omega, omega<0.01, 0.01)) %>%
         mutate(fitness=factor(fitness), 
                model=factor(model),
                generation=factor(generation))
df <- df %>% mutate(model=reorder(model, match(model, names(model_map))))

# instructions for plotting multiple plots in a PDF
# http://stackoverflow.com/a/20502085/1135316


tmp = df %>%
       group_by(generation) %>%  
      do( plot=plot.fitness.factors.log(.) + 
          ggtitle(sprintf("Generation %s", levels(.$generation)[unique(.$generation)])) +
            scale_shape_discrete(guide = FALSE) +
            theme(panel.grid.major = element_line(size = .5, color = "grey"),
                  axis.text.x = element_text(size=8,face="plain"),
                  axis.text.y = element_text(size=8,face="plain")
                  )
          )
plots <- tmp$plot
invisible(lapply(plots, print))








