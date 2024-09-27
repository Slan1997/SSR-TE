## This file contains codes to generate Figure 3 in the manuscript: 
# "Adaptive Sample Size Calculations for Rare Disease Clinical Trials Using a Totality of Evidence Approach"

## Figure 3: Empirical power associated with each scenario

## first, need to set path to the folder of result datasets.
setwd("~/Desktop/DSintern/manuscript material/Result datasets") 

library(pacman)
pacman::p_load(tidyr, dplyr,readr,stringr,ggplot2)

### combine OLS and permutation results
power1 = read_csv(paste0(getwd(),"/power_newdf_729.csv")) %>% mutate(
  tot_SS = n_control*(1+alloc_ratio),
  expected_SS = mean_new_N_ct*(1+alloc_ratio), # ESS
  max_SS = max_N_ct*(1+alloc_ratio) # MSS
)

power2 = read_csv(paste0(getwd(),"/power1e4_perm2e3_729.csv"))%>% mutate(
  tot_SS = n_control*(1+alloc_ratio),
  expected_SS = mean_new_N_ct*(1+alloc_ratio),
  max_SS = max_N_ct*(1+alloc_ratio)
)

power_long = bind_rows(power1%>%mutate(method='OLS'),
                       power2%>%mutate(method='Permutation') )
power_long

dt_plot =power_long %>% select(-scena) %>%
  mutate( n_setting = paste0("(",theta_k,",",alloc_ratio,")"), 
          `Target Power` = paste0(100*(1-beta),"%")
  )%>% 
  relocate(method,n_setting,`Target Power`,1) %>% #mutate(Power=1-beta) %>%
  rename(Method = method) %>%
  mutate(`Initially Planned SS` = factor(shrink_factor,levels=c(1,.6),labels=c("Original","Underestimated"))) %>%
  #relocate(Power,.before=beta) %>%
  mutate_at(c(1:6),as.factor)%>%
  mutate(Adaption = case_when( SSR == "no" ~ 1,
                               SSR == "yes" ~ 2) ) %>%
  mutate(Adaption = factor(Adaption,levels=1:2,
                           labels=c("No SSR","SSR")))
str(dt_plot)


ggplot(dt_plot ,
       aes(x=n_setting, y=power*100,color=Adaption)) + 
  # geom_line(aes(lty=Method#group=Adaption,
  #               ),
  #           lwd=.3) + 
  geom_point(aes(shape = Method) ,size=2.5) +
  geom_hline(aes(yintercept = 100*(1-beta)), linetype="dashed", 
             color = "darkgray", lwd=.3) +
  # geom_hline(yintercept=90, linetype="dashed", 
  #            color = "dark grey", lwd=.3) +
  scale_shape_manual(name = "Method", 
                     values = c(1, 3), 
                     labels = c("OLS","Permutation")) +
  scale_color_manual(name = "Adaption", 
                     values = c("black", "red"), 
                     labels = c("No SSR","SSR")) +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=10#,angle = 90, vjust = 0.5, hjust=1),
        ),
        panel.background = element_rect(fill = "white",
                                        #colour = "lightblue",
                                        linewidth = 0.2, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.2, linetype = 'solid',
                                        colour = "lightgrey") ,
        strip.text = element_text(size=10),
        strip.background = element_rect(fill="white",linetype = 'solid',colour = "black"),
        
        legend.position = c(.24,.15),
        legend.background = element_rect(fill="lightgrey", 
                                         linewidth=1, linetype="solid"),
        legend.box = "horizontal"
        # panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
        #                                 colour = "grey")
  ) +
  #scale_color_brewer(palette = 'R4') +
  #facet_wrap(~SSR,ncol=2, labeller = "label_both") +
  scale_y_continuous(breaks=seq(50, 100,5),limits = c(50,100)
  ) +
  labs(x = "(Expected Effect Size, Allocation Ratio)",#"Expect Cohen's d for Single Endpoint",
       y = "Empirical Power (%)") +
  facet_grid(`Target Power`~`Initially Planned SS`,labeller = label_both)

