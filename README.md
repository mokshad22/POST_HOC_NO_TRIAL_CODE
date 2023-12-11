# POST_HOC_NO_TRIAL_CODE
R code for analysis and plotting using Bayesian Framework


library(bayestestR)
library(brms)
library(tidybayes)
library(ggdist)
library(ggpubr)
library(sjPlot)
library(tableone)
library(tidycmprsk)
library(survival)
library(rms)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(survminer)
library (haven)
library(knitr)
library(kableExtra)
library(scales)
library(patchwork)

####Load datasets for WHITES and AA from the directory


###Model for WHITES ###
fit_pf_48_white <- brm(pf_48 ~  Group  + pf_baseline + Apache_II_baselinez + Age_baselinez +age2
                       + BMI_baselinez + Male_baseline, 
                       family = gaussian(),
                       prior = c(
                         prior(normal(0, 198), class = "b"), #2.5 x baseline  SD of pf
                         prior(normal(183.2,183.2), class = "Intercept"),
                         prior(student_t(3, 0, 77.2), class = "sigma")),seed = 123,
                       data = set_white)
summary(fit_pf_48_white)



### create two versions of our dataset, one with all Control and one with all NO ###
df_sim_t0 <- set_white %>% select(c(Age_baselinez,age2 , Apache_II_baselinez,Male_baseline,
                                    Race_baseline,BMI_baselinez,site_x,record_id)) |>
  mutate(Group = "Control") |> na.omit()

df_sim_t1 <- set_white %>% select(c(Age_baselinez,age2,Apache_II_baselinez,Male_baseline,
                                    Race_baseline,BMI_baselinez,site_x,record_id)) |>
  mutate(Group = "Nitric Oxide") |> na.omit()



pp_t0 <- posterior_epred(fit_pf_48_white, newdata = df_sim_t0)

pp_t1 <- posterior_epred(fit_pf_48_white, newdata = df_sim_t1)

####Main vector : Difference in pf_ratio between two Treatment arms in WHITES ##############
diff <- pp_t1 - pp_t0

###Calculating proportions of WHITES having a Treatment effct greater than 0
sor<- sort(diff)
sor_group_fit_pf_48_gt_0<- sor[sor >0]
num_pf_48<-length(sor_group_fit_pf_48_gt_0)
fit_pf_48_gt_0 <- num_pf_48 /4000
fit_pf_48_gt_0




fit_pf_48_AA <- brm(pf_48 ~  Group  + pf_baseline + Apache_II_baselinez + Age_baselinez +age2
                       + BMI_baselinez + Male_baseline, 
                       family = gaussian(),
                       prior = c(
                         prior(normal(0, 198), class = "b"), #2.5 x baseline  SD of pf
                         prior(normal(183.2,183.2), class = "Intercept"),
                         prior(student_t(3, 0, 77.2), class = "sigma")),seed = 123,
                       data = set_AA)
summary(fit_pf_48_AA)



# create two versions of our dataset, one with all Control and one with all NO for African American (AA)###
df_sim_t0_AA <- set_AA %>% select(c(Age_baselinez,age2 , Apache_II_baselinez,Male_baseline,
                                    Race_baseline,BMI_baselinez,site_x,record_id)) |>
  mutate(Group = "Control") |> na.omit()

df_sim_t1_AA <- set_AA %>% select(c(Age_baselinez,age2,Apache_II_baselinez,Male_baseline,
                                    Race_baseline,BMI_baselinez,site_x,record_id)) |>
  mutate(Group = "Nitric Oxide") |> na.omit()


pp_t0_AA <- posterior_epred(fit_pf_48_AA, newdata = df_sim_t0_AA)

pp_t1_AA <- posterior_epred(fit_pf_48_AA, newdata = df_sim_t1_AA)

diff_AA <- pp_t1_AA - pp_t0_AA



###Calculating proportions of AA having a Treatment effct greater than 0
sor_b<- sort(diff_AA)
sor_group_fit_pf_48_gt_0<- sor_b[sor_b >0]
num_pf_48<-length(sor_group_fit_pf_48_gt_0)
fit_pf_48_gt_0 <- num_pf_48 /4000
fit_pf_48_gt_0









############# PLOTTING CODE###############

my_col = "darkgray"
d<- density(sor)
MyDF_w <- as.data.frame(cbind(x=d$x,y=d$y))

summary(d$x)
summary(d$y)
shade_curve <- function(MyDF_w, zstart, zend, fill = "red", alpha = .5){
  geom_area(data = subset(MyDF_w, x >= 0
                          & x < 83.95),
            aes(y=y), fill = fill, color = NA, alpha = alpha)
}


e<- density(sor_b)
MyDF_b <- as.data.frame(cbind(x=e$x,y=e$y))
summary(e$x)
summary(e$y)
shade_curve_b <- function(MyDF_b, zstart, zend, fill = "red", alpha = .5){
  geom_area(data = subset(MyDF_b, x >= 0
                          & x < 237),
            aes(y=y), fill = fill, color = NA, alpha = alpha)
}




ggplot(data=MyDF_w,aes(x = x, y = y)) + geom_line(color="blue") +
  shade_curve(MyDF = MyDF_w, zstart = 0, zend = 76.52, fill = my_col, alpha = .5) +
  theme_classic() +scale_x_continuous(breaks=c(-100,-50, 0 , 50 ,100, 150,200,250))+
  ylab("Density") + xlab("Mean difference (Nitric Oxide - Control)") +
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "darkgreen", size=1)+ ggtitle("Change in P/F ratio after 48 hours Black White") + ggeasy::easy_center_title()+
  geom_line(data=MyDF_b,aes(x = x, y = y),color="red") +
  shade_curve_b(MyDF = MyDF_b, zstart = -108.21, zend = 236.68, fill = my_col, alpha = .5)
geom_area(data = subset(MyDF_b, x >= 0
                        & x < 237),
          aes(y=y), fill = my_col, color = NA, alpha = 0.5)






