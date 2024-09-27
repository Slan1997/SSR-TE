## This file contains codes to generate Table 1 in the manuscript: 
# "Adaptive Sample Size Calculations for Rare Disease Clinical Trials Using a Totality of Evidence Approach"

## Table 1: Summary of simulation results of power assessment across scenarios

## first, need to set path to the folder of result datasets.
setwd("~/Desktop/DSintern/manuscript material/Result datasets") 

library(pacman)
pacman::p_load(tidyr, dplyr,readr,stringr,ggplot2)


power1_wide_noSSR = power1 %>% filter(SSR=="no") %>% 
  select(-SSR) %>%
  rename(power_OLS_noSSR = power,
         expected_SS_OLS_noSSR =expected_SS ,
         max_SS_OLS_noSSR = max_SS ) 

power1_wide_SSR = power1 %>% filter(SSR=="yes") %>% 
  select(-SSR) %>%
  transmute(power_OLS_SSR = power,
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
         expected_SS_perm_noSSR =expected_SS ,
         max_SS_perm_noSSR = max_SS ) 

power2_wide_SSR = power2 %>% filter(SSR=="yes") %>% 
  select(-SSR) %>%
  transmute(power_perm_SSR = power,
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