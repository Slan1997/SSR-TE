## This file contains codes to generate Figure 2 in the manuscript: 
# "Adaptive Sample Size Calculations for Rare Disease Clinical Trials Using a Totality of Evidence Approach"

## Figure 2: Empirical type I error associated with each scenario

## first, need to set path to the folder of result datasets.
setwd("~/Desktop/DSintern/manuscript material/Result datasets") 

library(pacman)
pacman::p_load(tidyr, dplyr,readr,stringr,ggplot2)

# Empirical type I error
### combine OLS and permutation results
dt_plot0 = read_csv(paste0(getwd(),"/alpha_newdf729.csv")) %>% 
  filter(n_control<1200) %>% mutate(Method='OLS')

dt_plot1 = read_csv(paste0(getwd(),"/1e4debug_perm2e3_729.csv")) %>% 
  filter(n_control<1200)%>% mutate(Method='Permutation') %>% select(-scena)

dt_plot =bind_rows(dt_plot0,dt_plot1) %>%
  #mutate(n_setting = paste0("N_C=",n_control,",N_T=",alloc_ratio*n_control))%>% 
  mutate(n_setting = paste0("(",n_control*(1+alloc_ratio),",",alloc_ratio,")"))%>% 
  relocate(n_setting,Method,1) %>% #mutate(Power=1-beta) %>%
  #relocate(Power,.before=beta) %>%
  mutate_at(2:5,as.factor) %>%
  mutate(n_setting = factor(n_setting,
                            levels = paste0(rep("(",length(unique(dt_plot0$n_control))*2), 
                                            rep(sort(unique(dt_plot0$n_control)),each=2)*(1+rep(1:2,length(unique(dt_plot0$n_control)))),
                                            rep(",",length(unique(dt_plot0$n_control))*2),
                                            rep(1:2,length(unique(dt_plot0$n_control))),
                                            rep(")",length(unique(dt_plot0$n_control))*2)
                            ))
  )%>%
  mutate(Adaption = case_when( SSR == "no" ~ 1,
                               SSR == "yes" & beta==.1 ~ 2,
                               SSR == "yes" & beta==.2 ~ 3) ) %>%
  mutate(Adaption = factor(Adaption,levels=1:3,
                           labels=c("No SSR","SSR (90% Target Power)","SSR (80% Target Power)"
                           )
  ))
str(dt_plot)

######## theta_k==0
ggplot(dt_plot ,
       aes(x=n_setting, y=alpha*100,color=Adaption)) + 
  geom_line(aes(group=Adaption),
            lwd=.3) + 
  geom_point(aes(shape = Adaption) ) +
  geom_hline(yintercept=2.5, linetype="dashed", 
             color = "darkgrey", lwd=.3) +
  scale_shape_manual(name = "Adaption", 
                     values = c(1, 2,3), 
                     labels = c("No SSR","SSR (90% Target Power)","SSR (80% Target Power)"
                     )) +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1,size=9),
        axis.text.y = element_text(size=10),
        panel.background = element_rect(fill = "white",
                                        #colour = "lightblue",
                                        size = 0.2, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "lightgrey") ,
        strip.background = element_rect(fill="white",linetype = 'solid',colour = "black"),
        strip.text = element_text(size=10),
        legend.position = c(.2,.75),
        legend.background = element_rect(fill="lightgrey", 
                                         size=1, linetype="solid")
        # panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
        #                                 colour = "grey")
  ) +
  scale_color_brewer(palette = 'Set1') +
  #facet_wrap(~SSR,ncol=2, labeller = "label_both") +
  scale_y_continuous(breaks=seq(2.2, 3,0.1)#,limits = c(.023,.026)
  ) +
  labs(x = "(Initially Planned Sample Size, Allocation Ratio)",
       y = "Empirical Type I Error (%)") +
  facet_grid(~Method)
