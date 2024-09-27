## This file contains codes to generate Figure 1 in the manuscript: 
# "Adaptive Sample Size Calculations for Rare Disease Clinical Trials Using a Totality of Evidence Approach"

## Figure 1: Examples of one-sided efficacy boundaries using the 
# Lan-DeMets O'Brien-Fleming-type alpha spending function, with a 
# total number of stages ranging from 2 to 4 and a type I error level of $\alpha = 0.025$.

## first, need to set path to the folder of result datasets.
setwd("~/Desktop/DSintern/manuscript material/Result datasets") 

library(pacman)
pacman::p_load(tidyr, dplyr,readr,stringr,ggplot2)

############################################### Figure 1
# Examples of one-sided efficacy boundaries using the Lan-DeMets O'Brien-Fleming-type 
# alpha spending function

library(gsDesign)

## For 2-stage, what are the boundaries
alpha <- 0.025 # 1-sided Type I error
k <- 2 # Number of planned analyses
test.type <- 1 # 1=one-sided (we didn't consider futility bounds heree)
sfu <- sfLDOF # spending function by Lan and DeMets (Lan and DeMets 1983) 
x <- gsDesign(
  k = k, 
  test.type = test.type,
  alpha = alpha,
  timing =1/k,
  sfu = sfu
)
C1 = x$upper$bound[1]
C2 = x$upper$bound[2]
C1;C2

#### when the total number of stages ranges from 2 to 4 
alpha <- 0.025 # 1-sided Type I error
test.type <- 1 # Asymmetric bound design with non-binding futility bound
sfu <- sfLDOF # O'Brien-Fleming spending function for alpha-spending
bounds = list()
for (k in 2:4){
  timing <- 1/k # information fraction at interim analyses
  x <- gsDesign(
    k = k,
    test.type = test.type,
    alpha = alpha,
    timing = timing,
    sfu = sfu
  )
  bounds = c(bounds,list(x$upper$bound))
}
bounds

dt_bounds = tibble(x=c((1:2)/2,(1:3)/3,(1:4)/4),
                   bounds = do.call(c,bounds),
                   `Number of Stages` = factor(c(2,2,3,3,3,4,4,4,4),
                                               levels=2:4,
                                               labels=paste(2:4, "Stages")))

# plot the figure
ggplot(dt_bounds, aes(x=x, y=bounds,group=`Number of Stages`,color=`Number of Stages`)) +
  geom_point(aes(shape=`Number of Stages`),size=2) +
  geom_line(aes(lty =`Number of Stages` ),lwd=.6)+
  scale_linetype_manual(values=c("solid","twodash", "dotted")) +
  scale_shape_manual(values=1:3) +
  scale_color_brewer(palette = 'Dark2') +
  scale_x_continuous(breaks=c(1/4,1/3,1/2,2/3,3/4,1),
                     labels=c("1/4","1/3","1/2","2/3","3/4","1")#, limits = c(0,1)
  ) +
  theme_bw() + 
  theme(legend.position=c(.78,.78),
        legend.background = element_rect(fill="lightgrey", 
                                         linewidth=1, linetype="solid"),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.title.align=0.5) + 
  ylab("One-Sided Boundaries") + 
  xlab("Information Fraction")
