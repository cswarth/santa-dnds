


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
    plot.fitness.factors(df) + coord_trans(y="log10") + scale_y_continuous(breaks = breaks)
}


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
      facet_wrap(~ model, ncol=4, drop=TRUE, scales='free_y', labeller = labeller(model=model_map)) +
      scale_x_continuous(limits = c(0, 2)) +
      scale_colour_brewer(type='qual', palette="Paired", name="Fitness Strength", breaks=c(3, 2, 1),
                          labels=c("high", "middle", "low")) +
      scale_shape_discrete(name="Generation")+
      xlab(NULL)
    
    plot_tab <- ggplotGrob(my_plot)

    # elimintae the X-axis from all the facet plots.
    plot_filtered <- gtable_filter(plot_tab, "^(?!axis_b)",trim=FALSE, perl=TRUE)
   # grid.newpage()
  #  grid.draw(plot_filtered)
}

plot.fitness.curves <- function(df) {
  df <- df %>% 
    group_by(model) %>% 
    mutate(flevel=as.numeric(factor(fitness))) %>% 
    ungroup() %>% 
    mutate(flevel=factor(flevel))
  
  df %>% 
    mutate(flevel=as.numeric(levels(flevel))[flevel]) %>%
    ggplot(aes(flevel, omega, color=model, group=replicate)) +
    theme_bw() + 
    geom_line(size=1.4) +
    geom_point(size=1.6)
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
    geom_point(size=1.6) +
    ylim(0,2)
  
}

df <- read_csv('build/results.csv')
df <- df %>% filter(complete.cases(.)) %>%
        mutate(omega=replace(omega, omega<0.01, 0.01)) %>%
         mutate(fitness=factor(fitness), 
                model=factor(model),
                generation=factor(generation))

# instructions for plotting multiple plots in a PDF
# http://stackoverflow.com/a/20502085/1135316


tmp = df %>% group_by(generation) %>%  
      do( plot=plot.fitness.factors.log(.) + 
          ggtitle(sprintf("Generation %s", levels(.$generation)[unique(.$generation)])) +
            scale_shape_discrete(guide = FALSE) 
          )
plots <- tmp$plot
pdf("factors.pdf", width=7, height=5)
invisible(lapply(plots, print))
dev.off()

ggsave("factors.pdf", width=8.5, height=6, units="in", marrangeGrob(grobs = plots, nrow=1, ncol=1))

tmp <- df %>% group_by(generation) %>% do(plot=plot.fitness.facets(.))
plots <- tmp$plot
pdf("facets.pdf")
invisible(lapply(plots, print))
dev.off()







