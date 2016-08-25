


library(readr)
library(dplyr)
library(ggplot2)
library(gtable)
library(grid)

plot.fitness.facets <- function(df) {
    # a copy of Hadley's gtable_filter function adds a `perl` keyword arg 
    # to support negated regex.
    # https://github.com/hadley/gtable/blob/master/R/filter.r
    gtable_filter <- function(x, pattern, fixed = FALSE, trim = TRUE, perl=FALSE) {
      matches <- grepl(pattern, x$layout$name, fixed = fixed, perl=perl)
      x$layout <- x$layout[matches, , drop = FALSE]
      x$grobs <- x$grobs[matches]
      if (trim) x <- gtable_trim(x)
      x
    }
    
    model_map <- function(value) {
      c("noselection" = "Neutral",
        "empiricalvalues" = "EmpiricalValue",
        "purifyingchem" = "ChemicalAffinity",
        "homoresidue" = "Empirical\nHomoresidue")[value]
    }
    
    df <- df %>% mutate(one=1)
    df <- df %>% 
            group_by(model) %>% 
            mutate(flevel=as.numeric(factor(fitness))) %>% 
            ungroup() %>% 
            mutate(flevel=factor(flevel))
    my_plot <- ggplot(df, aes(one, omega, color=flevel, shape=generation)) +
      theme_bw() +
      geom_point(position = position_jitter(w = 0.25), size=3) +
      facet_wrap(~ model, drop=TRUE, scales='free_y', labeller = labeller(model=model_map)) +
      scale_x_continuous(limits = c(0, 2)) +
      scale_color_discrete(name="Fitness Strength", breaks=c(3, 2, 1),
                           labels=c("high", "middle", "low"))+
      scale_shape_discrete(name="Generation")+
      xlab(NULL)
    
    plot_tab <- ggplotGrob(my_plot)
    
    plot_filtered <- gtable_filter(plot_tab, "^(?!axis_b)",trim=FALSE, perl=TRUE)
    grid.newpage()
    grid.draw(plot_filtered)
}

plot.mean.fitness.curves<-function(df) {
    df <- df %>% 
    group_by(model) %>% 
    mutate(flevel=as.numeric(factor(fitness))) %>% 
    ungroup() %>% 
    mutate(flevel=factor(flevel))
  
    df %>% filter(generation==20000) %>% 
    group_by(model, flevel) %>% 
    summarize(omega = mean(omega)) %>% 
    ungroup() %>% mutate(flevel=as.numeric(levels(flevel))[flevel]) %>%
    ggplot(aes(flevel, omega, color=model)) +
    theme_bw() + 
    geom_line(size=1.4) +
    geom_point(size=1.6)
    
}

plot.fitness.curves<-function(df) {
  df <- df %>% 
    group_by(model) %>% 
    mutate(flevel=as.numeric(factor(fitness))) %>% 
    ungroup() %>% 
    mutate(flevel=factor(flevel))
  
  df %>% filter(generation==500) %>% 
    filter(model=='homoresidue') %>%
    mutate(flevel=as.numeric(levels(flevel))[flevel]) %>%
    ggplot(aes(flevel, omega, color=model, group=replicate)) +
    theme_bw() + 
    geom_line(size=1.4) +
    geom_point(size=1.6)
  
}

df <- read_csv('build/results.csv')
df <- df %>% filter(complete.cases(.)) %>%
         filter(omega < 10) %>%
         mutate(fitness=factor(fitness), 
                model=factor(model),
                generation=factor(generation))

plot.fitness.facets(df)

plot.mean.fitness.curves(df)


