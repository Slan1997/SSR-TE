## This file contains codes to generate figures and table in the manuscript: 
# "Adaptive Sample Size Calculations for Rare Disease Clinical Trials Using a Totality of Evidence Approach"

## first, need to set path to the folder of result datasets.
setwd("~/DSintern/manuscript material/Result datasets") 

library(pacman)
pacman::p_load(tidyr, dplyr,readr,stringr,ggplot2)

############################################### Figure 1
# Examples of one-sided efficacy boundaries using the Lan-DeMets O'Brien-Fleming-type 
# alpha spending function

library(gsDesign)

## For 2-stage, what are theboundaries
alpha <- 0.025 # 1-sided Type I error
k <- 2 # Number of planned analyses
test.type <- 1 # Asymmetric bound design with non-binding futility bound
sfu <- sfLDOF
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
                                         size=1, linetype="solid"),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.title.align=0.5) + 
  ylab("One-Sided Boundaries") + 
  xlab("Information Fraction")


############################################### Figure 2
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

############################################### Figure 3
# Empirical power

power1 = read_csv(paste0(getwd(),"/power_newdf_729.csv")) %>% mutate(
  tot_SS = n_control*(1+alloc_ratio),
  expected_SS = mean_new_N_ct*(1+alloc_ratio),
  max_SS = max_N_ct*(1+alloc_ratio)
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
  mutate( n_setting = paste0("(",theta_k,",",alloc_ratio,")"), # paste0("ES=",theta_k,",r=",alloc_ratio), #,",N=",tot_SS
          `Target Power` = paste0(100*(1-beta),"%")# paste0(expression(theta),"=",theta_k,",r=",alloc_ratio#,",N=",tot_SS
          #                  )
  )%>% 
  relocate(method,n_setting,`Target Power`,1) %>% #mutate(Power=1-beta) %>%
  rename(Method = method) %>%
  mutate(`Initially Planned SS` = factor(shrink_factor,levels=c(1,.6),labels=c("Original","Underestimated"))) %>%
  #relocate(Power,.before=beta) %>%
  mutate_at(c(1:6),as.factor)%>%
  # lv = sort(unique(dt_plot$n_setting))
  # #gsub("r=0.5","r=2",lv)
  # dt_plot  = dt_plot %>%
  #   mutate(n_setting = factor(n_setting,
  #                             levels = lv
  #                             )
  #   ) %>%
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
                                        size = 0.2, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                        colour = "lightgrey") ,
        strip.text = element_text(size=10),
        strip.background = element_rect(fill="white",linetype = 'solid',colour = "black"),
        
        legend.position = c(.24,.15),
        legend.background = element_rect(fill="lightgrey", 
                                         size=1, linetype="solid"),
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



############################################### Table 1
# Empirical power

power1_wide_noSSR = power1 %>% filter(SSR=="no") %>% 
  select(-SSR) %>%
  rename(power_OLS_noSSR = power,
         # mean_new_N_ct_OLS_noSSR = mean_new_N_ct,
         # max_new_N_ct_OLS_noSSR= max_N_ct ,
         #tot_SS_OLS =tot_SS,
         expected_SS_OLS_noSSR =expected_SS ,
         max_SS_OLS_noSSR = max_SS ) 

power1_wide_SSR = power1 %>% filter(SSR=="yes") %>% 
  select(-SSR) %>%
  transmute(power_OLS_SSR = power,
            # mean_new_N_ct_OLS_SSR = mean_new_N_ct,
            # max_new_N_ct_OLS_SSR= max_N_ct ,
            #tot_SS_OLS =tot_SS,
            expected_SS_OLS_SSR =expected_SS,
            max_SS_OLS_SSR = max_SS) 

OLS_wide = bind_cols(power1_wide_noSSR,power1_wide_SSR ) %>%
  arrange(desc(beta),desc(shrink_factor),theta_k,alloc_ratio) %>%
  transmute(beta,shrink_factor, theta_k,alloc_ratio,tot_SS,
            power_OLS_noSSR,expected_SS_OLS_noSSR=ceiling(expected_SS_OLS_noSSR),
            power_OLS_SSR,expected_SS_OLS_SSR=ceiling(expected_SS_OLS_SSR),max_SS_OLS_SSR)

power2_wide_noSSR = power2 %>% filter(SSR=="no") %>% 
  select(-SSR) %>%
  rename(power_perm_noSSR = power,
         # mean_new_N_ct_OLS_noSSR = mean_new_N_ct,
         # max_new_N_ct_OLS_noSSR= max_N_ct ,
         #tot_SS_OLS =tot_SS,
         expected_SS_perm_noSSR =expected_SS ,
         max_SS_perm_noSSR = max_SS ) 

power2_wide_SSR = power2 %>% filter(SSR=="yes") %>% 
  select(-SSR) %>%
  transmute(power_perm_SSR = power,
            # mean_new_N_ct_OLS_SSR = mean_new_N_ct,
            # max_new_N_ct_OLS_SSR= max_N_ct ,
            #tot_SS_OLS =tot_SS,
            expected_SS_perm_SSR =expected_SS,
            max_SS_perm_SSR = max_SS) 


Perm_wide = bind_cols(power2_wide_noSSR,power2_wide_SSR ) %>%
  arrange(desc(beta),desc(shrink_factor),theta_k,alloc_ratio) %>%
  transmute(beta,shrink_factor, theta_k,alloc_ratio,tot_SS,
            power_perm_noSSR,expected_SS_perm_noSSR=ceiling(expected_SS_perm_noSSR),
            power_perm_SSR,expected_SS_perm_SSR=ceiling(expected_SS_perm_SSR),max_SS_perm_SSR)

power_wide = full_join(OLS_wide,Perm_wide,
                       by=c("beta","shrink_factor","theta_k","alloc_ratio","tot_SS")) %>%
  transmute(`Target Power` = 1-beta,
            `IPSS Type` = factor(shrink_factor,labels = c("Underest.","Original")),
            bar_thata_star=theta_k,
            r = alloc_ratio,
            N_total = tot_SS,
            power_OLS_noSSR,expected_SS_OLS_noSSR,
            power_OLS_SSR,expected_SS_OLS_SSR,max_SS_OLS_SSR,
            power_perm_noSSR,expected_SS_perm_noSSR,
            power_perm_SSR,expected_SS_perm_SSR,max_SS_perm_SSR
  )

power_wide %>% write_csv("power_table.csv")

