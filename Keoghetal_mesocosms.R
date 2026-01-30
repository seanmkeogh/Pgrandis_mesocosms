library(lmerTest)
library(cowplot)
library(Hmisc)
library(viridis)
library(readxl)
library(car)
library(scales)
library(multcomp)
library(MuMIn)
library(performance)
library(MASS)
library(patchwork)
library(tidyverse)
library(grid)
library(boot)
library(lme4)

setwd('~/Box/Myprojects/Pgrandis_baskets/')
path_figs<-('figures/ggplot')

#################################
## YEAR 1 - DATASET ASSEMBLY
#################################
all_w_dead <- read_xlsx('Year1/data/grandis_baskets.xlsx',
                        sheet='all',na=c('NA','','No Data','Unknown'))
# Remove the five dead individuals
all<-all_w_dead %>% filter(is.na(mark) | mark == 'black' | mark == 'white')

# Change row number to distance from aerators
all<-all %>% mutate(dist_m = ifelse(row ==1, 0, ifelse(row == 2, 3.05, 
                                                       ifelse(row ==3, 9.144,
                                                              ifelse(row ==4, 33.528,
                                                                     ifelse(row==5, 99.06,195.072))))))

pre<- all %>% filter(time == 'pre')
post<- all %>% filter(time == 'post')

# Take averages for each basket
pre_avg<-pre %>% group_by(basket) %>% 
  dplyr::summarize(pre_length = mean(length, na.rm = T),
            pre_height = mean(height, na.rm = T),
            pre_width = mean(width, na.rm = T))
# Not calculate growth (change) in measurements by basket averages to each individual
yr1_df <- pre_avg %>% left_join(post %>% select(-time,-mark), by='basket')%>% 
  mutate(length_growth = length - pre_length,
         height_growth = height - pre_height,
         width_growth = width - pre_width)

# Make sure categorical variables are categorical
yr1_df$row <- as.factor(yr1_df$row)
yr1_df$basket <- as.factor(yr1_df$basket)
# Change sand substrate so models correspond to substrate size
yr1_df <- yr1_df %>%
  mutate(substrate = recode(substrate, "sand" = "asand"))
yr1_df$substrate <- as.factor(yr1_df$substrate)
# Change density so models correspond to increasing density
yr1_df <- yr1_df %>%
  mutate(density = recode(density, "low" = "alow"))
yr1_df$density <- as.factor(yr1_df$density)
# Change distance to Flow/aeration to Proximity to Flow
yr1_df$dist_m <- -1*yr1_df$dist_m

#################################
## Year 2 - DATASET ASSEMBLY
#################################
baskets<-read_xlsx('Year2/data/Year2_data.xlsx',
                   sheet='baskets',na=c('NA',''),
                   col_types = c('text','text','text','text','numeric','text','numeric','numeric','numeric','numeric','numeric',
                                 'numeric','numeric','numeric','numeric','numeric',
                                 'numeric','numeric','numeric','numeric','numeric','numeric','numeric','numeric'))
mussels<-read_xlsx('Year2/data/Year2_data.xlsx',
                   sheet='mussels',na=c('NA',''), col_types = c('text','numeric','numeric','numeric','numeric','text','text','text','numeric')
                   ) %>% filter(Condition == 'alive')

# assemble data
baskets <- baskets %>% dplyr::select (Basket, Group, Substrate, Vert_pos, Depth, Location, avg_Temp, avg_DO, avg_pH, Flow)

df<-right_join(baskets,mussels, by = 'Basket')

may <- df %>% filter(Time == 'May') %>% 
  dplyr::select(Label,Length,Height,Width,Mass,Shell_mass_g)
may <- may %>% mutate(
  pre_Length = Length,
  pre_Height = Height,
  pre_Width = Width,
  pre_Mass = Mass) %>% dplyr::select(Label,pre_Length,pre_Height,pre_Width,pre_Mass)
oct <- df %>% filter(Time == 'October')

yr2_df<-oct %>% left_join(may,by = 'Label') %>% 
  mutate(Length_growth = Length - pre_Length,
         Height_growth = Height - pre_Height,
         Width_growth = Width - pre_Width,
         Mass_growth = Mass - pre_Mass)

# Make sure categorical variables are categorical
yr2_df$Basket <- as.factor(yr2_df$Basket)
yr2_df$Substrate <- as.factor(yr2_df$Substrate)
yr2_df$Group <- as.factor(yr2_df$Group)
yr2_df$Vert_pos <- as.factor(yr2_df$Vert_pos)

# Shell Thickness Index (STI) from Freeman and Byers 2006
yr2_df <- yr2_df %>% mutate(STI = 1000 * (Shell_mass_g/(Length * (((Height^2) + Width^2)^0.5) * pi/2)))

# STI taken from wet mass
yr2_df <- yr2_df %>% mutate(STI_wet = 1000 * (Mass/(Length * (((Height^2) + Width^2)^0.5) * pi/2)))

# Initial wet STI
yr2_df <- yr2_df %>% mutate(pre_STI_wet = 1000 * (pre_Mass/(pre_Length * (((pre_Height^2) + pre_Width^2)^0.5) * pi/2)))

# STI change
yr2_df <- yr2_df %>% mutate(STI_growth = STI_wet-pre_STI_wet)

# Soft tissue
yr2_df <- yr2_df %>% mutate(Soft_tissue = Mass-Shell_mass_g)

# Convert temp to C
yr2_df <- yr2_df %>% mutate(avg_Temp_c = (avg_Temp - 32) * (5/9))

# Convert Flow to m/s
yr2_df <- yr2_df %>% mutate(Flow_m = Flow * 0.3048)


# Save datasets
# saveRDS(yr1_df,file='data/yr1_df.data')
# saveRDS(yr2_df,file='data/yr2_df.data')
yr1_df <- readRDS(file='data/yr1_df.data')
yr2_df <- readRDS(file='data/yr2_df.data')




###############################################
#
#################################
#
### LINEAR MIXED EFFECTS MODELING ###
#
#################################
#
###############################################


#################################
#
### START: LOAD MODELS & BOOTSTRAPPING RESULTS
#
#################################

yr1.len <- readRDS(file='data/yr1.len.data')
yr1.len_boot_fixef_lmer <- readRDS(file='data/yr1.len_boot_fixef_lmer.data')
pred_sub_growth_BOOT <- readRDS(file='data/pred_sub_growth_BOOT.data')
pred_den_growth_BOOT <- readRDS(file='data/pred_den_growth_BOOT.data')

yr2.len <- readRDS(file='data/yr2.len.data')
yr2.len_boot_fixef_lmer <- readRDS(file='data/yr2.len_boot_fixef_lmer.data')
pred_flow_growth_BOOT <- readRDS(file='data/pred_flow_growth_BOOT.data')
pred_temp_growth_BOOT <- readRDS(file='data/pred_temp_growth_BOOT.data')

yr1.ht <- readRDS(file='data/yr1.ht.data')
yr1.ht_boot_fixef_lmer <- readRDS(file='data/yr1.ht_boot_fixef_lmer.data')
pred_den_ht_BOOT <- readRDS(file='data/pred_den_ht_BOOT.data')

yr2.ht <- readRDS(file='data/yr2.ht.data')
yr2.ht_boot_fixef_lmer <- readRDS(file='data/yr2.ht_boot_fixef_lmer.data')
pred_flow_ht_BOOT <- readRDS(file='data/pred_flow_ht_BOOT.data')
pred_temp_ht_BOOT <- readRDS(file='data/pred_temp_ht_BOOT.data')

yr1.wdt <- readRDS(file='data/yr1.wdt.data')
yr1.wdt_boot_fixef_lmer <- readRDS(file='data/yr1.wdt_boot_fixef_lmer.data')
pred_sub_wt_BOOT <- readRDS(file='data/pred_sub_wt_BOOT.data')
pred_den_wt_BOOT <- readRDS(file='data/pred_den_wt_BOOT.data')

yr2.wdt <- readRDS(file='data/yr2.wdt.data')
yr2.wdt_boot_fixef_lmer <- readRDS(file='data/yr2.wdt_boot_fixef_lmer.data')

yr2.sti <- readRDS(file='data/yr2.sti.data')
yr2.sti_boot_fixef_lmer <- readRDS(file='data/yr2.sti_boot_fixef_lmer.data')
pred_flow_sti_BOOT <- readRDS(file='data/pred_flow_sti_BOOT.data')
pred_temp_sti_BOOT <- readRDS(file='data/pred_temp_sti_BOOT.data')

#################################
#
### END: LOAD MODELS & BOOTSTRAPPING RESULTS
#
#################################



###################################
# SHELL LENGTH - GROWTH RATE
###################################

# Function for bootstrapping
fixef_fun <- function(fit) fixef(fit)
set.seed(17)

# YEAR 1 - GROWTH
# yr1.len <- lmer(length_growth ~ substrate + dist_m + density + pre_length + (1|basket), data=yr1_df, REML = T)
# saveRDS(yr1.len,file='data/yr1.len.data')
check_model(yr1.len, check = c("linearity", "homogeneity", "qq", "normality"))
summary(yr1.len)
vc <- as.data.frame(VarCorr(yr1.len))
var_basket <- vc$vcov[vc$grp == "basket"]
var_resid <- vc$vcov[vc$grp == "Residual"]
ICC <- var_basket / (var_basket + var_resid)
# yr1.len_boot_fixef_lmer <- bootMer(yr1.len, fixef_fun,nsim = 10000, use.u = FALSE, type = "parametric")
# saveRDS(yr1.len_boot_fixef_lmer,file='data/yr1.len_boot_fixef_lmer.data')
yr1.len_ci_lmer <- apply(
  yr1.len_boot_fixef_lmer$t,
  2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
yr1.len_ci_lmer <- t(yr1.len_ci_lmer)
colnames(yr1.len_ci_lmer) <- c("CI_lower", "CI_upper")
data.frame(
  Estimate = fixef(yr1.len),
  CI_lower = yr1.len_ci_lmer[, "CI_lower"],
  CI_upper = yr1.len_ci_lmer[, "CI_upper"])
# ^ substrate & density non-zero 95% CI!


# # YEAR 1 - GROWTH - SUBSTRATE PLOTTING
# New data - keeping other fixed effects constant
pred_sub_growth_data <- expand.grid(
  substrate = levels(factor(yr1_df$substrate)),
  dist_m = mean(yr1_df$dist_m, na.rm = TRUE),
  density = levels(factor(yr1_df$density))[1],
  pre_length = mean(yr1_df$pre_length, na.rm = TRUE))
# Predictions
pred_sub_growth_data$fit <- predict(yr1.len, newdata = pred_sub_growth_data, re.form = NA)
# Bootstrap predictions
pred_sub_growth_FUN <- function(fit) {
  predict(fit, newdata = pred_sub_growth_data, re.form = NA)
}
# Run parametric bootstrapping
# pred_sub_growth_BOOT <- bootMer(yr1.len, pred_sub_growth_FUN,
#   nsim = 10000, use.u = FALSE, type = "parametric")
# saveRDS(pred_sub_growth_BOOT,file='data/pred_sub_growth_BOOT.data')
# Lower bound of 95% CI
pred_sub_growth_data$lower <- apply(pred_sub_growth_BOOT$t, 2,
  quantile, probs = 0.025, na.rm = TRUE)
# Upper bound of 95% CI
pred_sub_growth_data$upper <- apply(pred_sub_growth_BOOT$t, 2,
  quantile, probs = 0.975, na.rm = TRUE)
boot_df <- data.frame(value = as.vector(pred_sub_growth_BOOT$t),
  substrate = rep(pred_sub_growth_data$substrate, each = nrow(pred_sub_growth_BOOT$t)))
yr1_sub_grow_plot <- ggplot() +
  geom_violin(data = boot_df, 
              aes(x = substrate, y = value, fill = substrate), trim = FALSE) +
  scale_fill_manual(values = c('#AF804F','#7F5112')) +
  geom_point(data = pred_sub_growth_data,
              aes(x = substrate, y = fit), size = 3) +
  geom_errorbar(data = pred_sub_growth_data,
              aes(x = substrate, ymin = lower, ymax = upper), width = 0.15) +
  scale_x_discrete(labels= c('Sand', 'Gravel')) +
  ylim(c(23,50.35)) +
  labs(x = "Substrate size", y = "Predicted \u0394 shell length (mm)") + 
  guides(fill="none") +
  theme_linedraw() +
  theme(plot.margin = margin(2, 2, 2, 2))
# Get Histo for plot above
sub_boot <- yr1.len_boot_fixef_lmer$t[, "substrategravel"]
sub_ci <- quantile(sub_boot, probs = c(0.025, 0.975), na.rm = TRUE)
histo <- ggplot(data.frame(sub_boot), aes(x = sub_boot)) +
  geom_density(fill = "grey80") +
  # xlim(c(-20,2)) +
  # ylim(c(0,.2)) +
  scale_y_continuous(limits=c(0, .2), breaks = seq(0, .2, by = .1)) +
  scale_x_continuous(limits=c(-20, 2), breaks = c(-20, -10, 0)) +
  geom_vline(xintercept = 0, linewidth = 0.7,
             linetype = "solid",color= "red") +
  geom_vline(xintercept = sub_ci,linewidth = 0.5,
             linetype = "solid",color= "blue",alpha=0.65) +
  labs(x = "Bootstrap estimates of Substrate effect",y = "Density") +
  theme_linedraw()
histo <- histo +
  theme_linedraw(base_size = 8) +
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.margin = margin(2, 2, 2, 2),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
# Combined plot for EXP I Substrate effect on Growth
yr1_sub_grow_plot_inset <- yr1_sub_grow_plot +
  inset_element(
    histo,
    left = 0.7,    # x-position of left edge (0–1)
    bottom = 0.7,  # y-position of bottom edge (0–1)
    right = 0.98,   # x-position of right edge
    top = 0.98,     # y-position of top edge
    align_to = "panel")

# # YEAR 1 - GROWTH - DENSITY PLOTTING
# New data - keeping other fixed effects constant
pred_den_growth_data <- expand.grid(
  substrate = levels(factor(yr1_df$substrate))[1],
  dist_m = mean(yr1_df$dist_m, na.rm = TRUE),
  density = levels(factor(yr1_df$density)),
  pre_length = mean(yr1_df$pre_length, na.rm = TRUE))
# Predictions
pred_den_growth_data$fit <- predict(yr1.len, newdata = pred_den_growth_data, re.form = NA)
# Bootstrap predictions
pred_den_growth_FUN <- function(fit) {
  predict(fit, newdata = pred_den_growth_data, re.form = NA)
}
# Run parametric bootstrapping
# pred_den_growth_BOOT <- bootMer(yr1.len, pred_den_growth_FUN,
#                                 nsim = 10000, use.u = FALSE, type = "parametric")
# saveRDS(pred_den_growth_BOOT,file='data/pred_den_growth_BOOT.data')
# Lower bound of 95% CI
pred_den_growth_data$lower <- apply(pred_den_growth_BOOT$t, 2,
                                    quantile, probs = 0.025, na.rm = TRUE)
# Upper bound of 95% CI
pred_den_growth_data$upper <- apply(pred_den_growth_BOOT$t, 2,
                                    quantile, probs = 0.975, na.rm = TRUE)
boot_df <- data.frame(value = as.vector(pred_den_growth_BOOT$t),
                      density = rep(pred_den_growth_data$density,
                                      each = nrow(pred_den_growth_BOOT$t)))
yr1_den_grow_plot <- ggplot() +
  geom_violin(data = boot_df, 
              aes(x = density, y = value, fill = density), trim = FALSE) +
  scale_fill_manual(values = c('gray90','gray40')) +
  geom_point(data = pred_den_growth_data,
             aes(x = density, y = fit), size = 3) +
  geom_errorbar(data = pred_den_growth_data,
                aes(x = density, ymin = lower, ymax = upper), width = 0.15) +
  scale_x_discrete(labels= c('Low - 10 indivs.', 'High - 50 indivs.')) +
  ylim(c(23,50.35)) +
  labs(x = "Mussel density", y = "Predicted \u0394 shell length (mm)") + 
  guides(fill="none") +
  theme_linedraw() +
  theme(plot.margin = margin(2, 2, 2, 2))
# Get Histo for plot above
den_boot <- yr1.len_boot_fixef_lmer$t[, "densityhigh"]
den_ci <- quantile(den_boot, probs = c(0.025, 0.975), na.rm = TRUE)
histo <- ggplot(data.frame(den_boot), aes(x = den_boot)) +
  geom_density(fill = "grey80") +
  # xlim(c(-20,2)) +
  # ylim(c(0,.2)) +
  scale_y_continuous(limits=c(0, .2), breaks = seq(0, .2, by = .1)) +
  scale_x_continuous(limits=c(-20, 2), breaks = c(-20, -10, 0)) +
  geom_vline(xintercept = 0, linewidth = 0.7,
             linetype = "solid",color= "red") +
  geom_vline(xintercept = den_ci,linewidth = 0.5,
             linetype = "solid",color= "blue",alpha=0.65) +
  labs(x = "Bootstrap estimates of Density effect",y = "Density") +
  theme_linedraw()
histo <- histo +
  theme_linedraw(base_size = 8) +
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.margin = margin(2, 2, 2, 2),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
# Combined plot for EXP I Substrate effect on Growth
yr1_den_grow_plot_inset <- yr1_den_grow_plot +
  inset_element(
    histo,
    left = 0.7,    # x-position of left edge (0–1)
    bottom = 0.7,  # y-position of bottom edge (0–1)
    right = 0.98,   # x-position of right edge
    top = 0.98,     # y-position of top edge
    align_to = "panel")


# YEAR 2 - GROWTH
# yr2.len <- lmer(Length_growth ~  Substrate + Flow_m + avg_Temp_c + pre_Length + (1|Basket), data=yr2_df, REML = T)
# saveRDS(yr2.len,file='data/yr2.len.data')
check_model(yr2.len, check = c("linearity", "homogeneity", "qq", "normality"))
summary(yr2.len)
vc <- as.data.frame(VarCorr(yr2.len))
var_basket <- vc$vcov[vc$grp == "Basket"]
var_resid <- vc$vcov[vc$grp == "Residual"]
var_basket / (var_basket + var_resid)
# yr2.len_boot_fixef_lmer <- bootMer(yr2.len, fixef_fun,nsim = 10000, use.u = FALSE, type = "parametric")
# saveRDS(yr2.len_boot_fixef_lmer,file='data/yr2.len_boot_fixef_lmer.data')
yr2.len_ci_lmer <- apply(
  yr2.len_boot_fixef_lmer$t,
  2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
yr2.len_ci_lmer <- t(yr2.len_ci_lmer)
colnames(yr2.len_ci_lmer) <- c("CI_lower", "CI_upper")
data.frame(
  Estimate = fixef(yr2.len),
  CI_lower = yr2.len_ci_lmer[, "CI_lower"],
  CI_upper = yr2.len_ci_lmer[, "CI_upper"])
# Flow and Temperature have non-zero effect on growth!

# # YEAR 2 - GROWTH - FLOW PLOTTING
# New data - keeping other fixed effects constant
flow_range <- seq(
  min(yr2_df$Flow_m, na.rm = TRUE),
  max(yr2_df$Flow_m, na.rm = TRUE),
  length.out = 100)
pred_flow_growth_data <- expand.grid(
  Flow_m = flow_range,
  Substrate = levels(factor(yr2_df$Substrate))[1],
  pre_Length = mean(yr2_df$pre_Length, na.rm = TRUE),
  avg_Temp_c = mean(yr2_df$avg_Temp_c, na.rm = TRUE))
# Predictions
pred_flow_growth_data$fit <- predict(yr2.len, newdata = pred_flow_growth_data, re.form = NA)
# Bootstrap predictions
pred_flow_growth_FUN <- function(fit) {
  predict(fit, newdata = pred_flow_growth_data, re.form = NA)
}
# Run parametric bootstrapping
# pred_flow_growth_BOOT <- bootMer(yr2.len, pred_flow_growth_FUN,
#                                 nsim = 10000, use.u = FALSE, type = "parametric")
# saveRDS(pred_flow_growth_BOOT,file='data/pred_flow_growth_BOOT.data')
# Lower bound of 95% CI
pred_flow_growth_data$lower <- apply(pred_flow_growth_BOOT$t, 2,
                                    quantile, probs = 0.025, na.rm = TRUE)
# Upper bound of 95% CI
pred_flow_growth_data$upper <- apply(pred_flow_growth_BOOT$t, 2,
                                    quantile, probs = 0.975, na.rm = TRUE)
boot_df <- data.frame(value = as.vector(pred_flow_growth_BOOT$t),
                      Flow_m = rep(pred_flow_growth_data$Flow_m,
                                    each = nrow(pred_flow_growth_BOOT$t)),
                      bootstrap_sample = paste0('bootstrap_',
                                                row_number(pred_flow_growth_BOOT$t)))
yr1_flow_grow_plot <- ggplot() +
  geom_line(data = boot_df, 
            aes(x = Flow_m, y = value, 
                group = bootstrap_sample),
            alpha = 0.05, size = 0.3, color='royalblue') +
  geom_line(data = pred_flow_growth_data,
            aes(x = Flow_m, y = fit),
            color = 'black',
            size = 1.5) +
  scale_x_continuous(limits = range(boot_df$Flow_m),
    expand = c(0, 0)) +
  ylim(c(18.9,50)) +
  labs(x = "Flow velocity (m/s)",
       y = "Predicted \u0394 shell length (mm)") +
  theme_linedraw() +
  theme(plot.margin = margin(2, 2, 2, 2))
max(pred_flow_growth_data$fit) - min(pred_flow_growth_data$fit)
# Get Histo for plot above
flow_boot <- yr2.len_boot_fixef_lmer$t[, "Flow_m"]
flow_ci <- quantile(flow_boot, probs = c(0.025, 0.975), na.rm = TRUE)
histo <- ggplot(data.frame(flow_boot), aes(x = flow_boot)) +
  geom_density(fill = "grey80") +
  # xlim(c(-2,54)) +
  # ylim(c(0,.3)) +
  scale_y_continuous(limits=c(0, .33), breaks = seq(0, .3, by = .15)) +
  scale_x_continuous(limits=c(-2, 54), breaks = c(0, 26, 54)) +
  geom_vline(xintercept = 0, linewidth = 0.7,
             linetype = "solid",color= "red") +
  geom_vline(xintercept = flow_ci,linewidth = 0.5,
             linetype = "solid",color= "blue",alpha=0.65) +
  labs(x = "Bootstrap estimates of Density effect",y = "Density") +
  theme_linedraw()
histo <- histo +
  theme_linedraw(base_size = 8) +
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.margin = margin(2, 2, 2, 2),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
# Combined plot for EXP II FLOW effect on Growth
yr1_flow_grow_plot_inset <- yr1_flow_grow_plot +
  inset_element(
    histo,
    left = 0.02,    # x-position of left edge (0–1)
    bottom = 0.7,  # y-position of bottom edge (0–1)
    right = 0.3,   # x-position of right edge
    top = 0.98,     # y-position of top edge
    align_to = "panel")

# # YEAR 2 - GROWTH - TEMP PLOTTING
# New data - keeping other fixed effects constant
temp_range <- seq(
  min(yr2_df$avg_Temp_c, na.rm = TRUE),
  max(yr2_df$avg_Temp_c, na.rm = TRUE),
  length.out = 100)
pred_temp_growth_data <- expand.grid(
  avg_Temp_c = temp_range,
  Substrate = levels(factor(yr2_df$Substrate))[1],
  pre_Length = mean(yr2_df$pre_Length, na.rm = TRUE),
  Flow_m = mean(yr2_df$Flow_m, na.rm = TRUE))
# Predictions
pred_temp_growth_data$fit <- predict(yr2.len, newdata = pred_temp_growth_data, re.form = NA)
# Bootstrap predictions
pred_temp_growth_FUN <- function(fit) {
  predict(fit, newdata = pred_temp_growth_data, re.form = NA)
}
# Run parametric bootstrapping
# pred_temp_growth_BOOT <- bootMer(yr2.len, pred_temp_growth_FUN,
#                                  nsim = 10000, use.u = FALSE, type = "parametric")
# saveRDS(pred_temp_growth_BOOT,file='data/pred_temp_growth_BOOT.data')
# Lower bound of 95% CI
pred_temp_growth_data$lower <- apply(pred_temp_growth_BOOT$t, 2,
                                     quantile, probs = 0.025, na.rm = TRUE)
# Upper bound of 95% CI
pred_temp_growth_data$upper <- apply(pred_temp_growth_BOOT$t, 2,
                                     quantile, probs = 0.975, na.rm = TRUE)
boot_df <- data.frame(value = as.vector(pred_temp_growth_BOOT$t),
                      avg_Temp_c = rep(pred_temp_growth_data$avg_Temp_c,
                                   each = nrow(pred_temp_growth_BOOT$t)),
                      bootstrap_sample = paste0('bootstrap_',
                                                row_number(pred_temp_growth_BOOT$t)))
yr1_temp_grow_plot <- ggplot() +
  geom_line(data = boot_df, 
            aes(x = avg_Temp_c, y = value, 
                group = bootstrap_sample),
            alpha = 0.05, size = 0.3, color='#DC143C') +
  geom_line(data = pred_temp_growth_data,
            aes(x = avg_Temp_c, y = fit),
            color = 'black',
            size = 1.5) +
  scale_x_continuous(limits = range(boot_df$avg_Temp_c),
                     expand = c(0, 0)) +
  ylim(c(18.9,50)) +
  labs(x = "Water temperature (°C)",
       y = "Predicted \u0394 shell length (mm)") +
  theme_linedraw() +
  theme(plot.margin = margin(2, 2, 2, 2))
max(pred_temp_growth_data$fit) - min(pred_temp_growth_data$fit)
# Get Histo for plot above
temp_boot <- yr2.len_boot_fixef_lmer$t[, "avg_Temp_c"]
temp_ci <- quantile(temp_boot, probs = c(0.025, 0.975), na.rm = TRUE)
histo <- ggplot(data.frame(temp_boot), aes(x = temp_boot)) +
  geom_density(fill = "grey80") +
  # xlim(c(-2,54)) +
  # ylim(c(0,.3)) +
  scale_y_continuous(limits=c(0, .33), breaks = seq(0, .3, by = .15)) +
  scale_x_continuous(limits=c(-2, 54), breaks = c(0, 26, 54)) +
  geom_vline(xintercept = 0, linewidth = 0.7,
             linetype = "solid",color= "red") +
  geom_vline(xintercept = temp_ci,linewidth = 0.5,
             linetype = "solid",color= "blue",alpha=0.65) +
  labs(x = "Bootstrap estimates of Temp effect",y = "Density") +
  theme_linedraw()
histo <- histo +
  theme_linedraw(base_size = 8) +
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.margin = margin(2, 2, 2, 2),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
# Combined plot for EXP II TEMP effect on Growth
yr1_temp_grow_plot_inset <- yr1_temp_grow_plot +
  inset_element(
    histo,
    left = 0.02,    # x-position of left edge (0–1)
    bottom = 0.7,  # y-position of bottom edge (0–1)
    right = 0.3,   # x-position of right edge
    top = 0.98,     # y-position of top edge
    align_to = "panel")

###
## PLOT SHELL GROWTH EXP I & II TOGETHER
###

left_panel <- (yr1_sub_grow_plot_inset / plot_spacer() /
                 yr1_den_grow_plot_inset) + 
  plot_layout(heights = c(1,0.05,1))
right_panel <- (yr1_flow_grow_plot_inset / plot_spacer() /
                  yr1_temp_grow_plot_inset) + 
  plot_layout(heights = c(1,0.05,1))
exp1_title <- textGrob(
  "Experiment I",
  gp = gpar(family = 'serif',fontface = "italic", fontsize = 11))
exp2_title <- textGrob(
  "Experiment II",
  gp = gpar(family = 'serif',fontface = "italic", fontsize = 11))
growth_figure <- ((wrap_elements(exp1_title) |
                     wrap_elements(exp2_title)) /
    ((left_panel | right_panel))) +
  plot_layout(
    heights = c(0.05, 1), widths = c(1,0.1,1)) 

growth_figure <- growth_figure &
  theme(text = element_text(family = "serif"))

ggsave(growth_figure,filename='growth_figure.png',path = paste0(path_figs),
       width = 6, height = 5.5, device='png', dpi=700)



###################################
# END OF SHELL LENGTH - GROWTH RATE
###################################



###################################
# SHELL HEIGHT - ROUNDNESS
###################################

# YEAR 1 - ROUNDNESS
# yr1.ht <- lmer(height ~ substrate + dist_m  + density + length + pre_height + pre_length + (1|basket), data=yr1_df, REML = T)
# saveRDS(yr1.ht,file='data/yr1.ht.data')
check_model(yr1.ht, check = c("linearity", "homogeneity", "qq", "normality"))
summary(yr1.ht)
vc <- as.data.frame(VarCorr(yr1.ht))
var_basket <- vc$vcov[vc$grp == "basket"]
var_resid <- vc$vcov[vc$grp == "Residual"]
var_basket / (var_basket + var_resid)
# yr1.ht_boot_fixef_lmer <- bootMer(yr1.ht, fixef_fun,nsim = 10000, use.u = FALSE, type = "parametric")
# saveRDS(yr1.ht_boot_fixef_lmer,file='data/yr1.ht_boot_fixef_lmer.data')
yr1.ht_ci_lmer <- apply(
  yr1.ht_boot_fixef_lmer$t,
  2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
yr1.ht_ci_lmer <- t(yr1.ht_ci_lmer)
colnames(yr1.ht_ci_lmer) <- c("CI_lower", "CI_upper")
data.frame(
  Estimate = fixef(yr1.ht),
  CI_lower = yr1.ht_ci_lmer[, "CI_lower"],
  CI_upper = yr1.ht_ci_lmer[, "CI_upper"])
# Density has creditable effect on shell height (roundness)!

# # YEAR 1 - ROUNDNESS - DENSITY PLOTTING
# New data - keeping other fixed effects constant
pred_den_ht_data <- expand.grid(
  substrate = levels(factor(yr1_df$substrate))[1],
  dist_m = mean(yr1_df$dist_m, na.rm = TRUE),
  density = levels(factor(yr1_df$density)),
  pre_length = mean(yr1_df$pre_length, na.rm = TRUE),
  pre_height = mean(yr1_df$pre_height, na.rm = TRUE),
  length = mean(yr1_df$length, na.rm = TRUE))
# Predictions
pred_den_ht_data$fit <- predict(yr1.ht, newdata = pred_den_ht_data, re.form = NA)
# Bootstrap predictions
pred_den_ht_FUN <- function(fit) {
  predict(fit, newdata = pred_den_ht_data, re.form = NA)
}
# Run parametric bootstrapping
# pred_den_ht_BOOT <- bootMer(yr1.ht, pred_den_ht_FUN,
#                                 nsim = 10000, use.u = FALSE, type = "parametric")
# saveRDS(pred_den_ht_BOOT,file='data/pred_den_ht_BOOT.data')
# Lower bound of 95% CI
pred_den_ht_data$lower <- apply(pred_den_ht_BOOT$t, 2,
                                    quantile, probs = 0.025, na.rm = TRUE)
# Upper bound of 95% CI
pred_den_ht_data$upper <- apply(pred_den_ht_BOOT$t, 2,
                                    quantile, probs = 0.975, na.rm = TRUE)
boot_df <- data.frame(value = as.vector(pred_den_ht_BOOT$t),
                      density = rep(pred_den_ht_data$density,
                                    each = nrow(pred_den_ht_BOOT$t)))
yr1_den_ht_plot <- ggplot() +
  geom_violin(data = boot_df, 
              aes(x = density, y = value, fill = density), trim = FALSE) +
  scale_fill_manual(values = c('gray90','gray40')) +
  geom_point(data = pred_den_ht_data,
             aes(x = density, y = fit), size = 3) +
  geom_errorbar(data = pred_den_ht_data,
                aes(x = density, ymin = lower, ymax = upper), width = 0.15) +
  scale_x_discrete(labels= c('Low - 10 indivs.', 'High - 50 indivs.')) +
  #ylim(c(23,50.35)) +
  labs(x = "Mussel density", y = "Predicted shell height (mm)") + 
  guides(fill="none") +
  theme_linedraw() +
  theme(plot.margin = margin(2, 2, 2, 2))
# Get Histo for plot above
den_boot <- yr1.ht_boot_fixef_lmer$t[, "densityhigh"]
den_ci <- quantile(den_boot, probs = c(0.025, 0.975), na.rm = TRUE)
histo <- ggplot(data.frame(den_boot), aes(x = den_boot)) +
  geom_density(fill = "grey80") +
  # xlim(c(-20,2)) +
  # ylim(c(0,.2)) +
  scale_y_continuous(limits=c(0, 1.6), breaks = seq(0, 1.6, by = .8)) +
  scale_x_continuous(limits=c(-2.2, 0.1), breaks = c(-2, -1, 0)) +
  geom_vline(xintercept = 0, linewidth = 0.7,
             linetype = "solid",color= "red") +
  geom_vline(xintercept = den_ci,linewidth = 0.5,
             linetype = "solid",color= "blue",alpha=0.65) +
  labs(x = "Bootstrap estimates of Density effect",y = "Density") +
  theme_linedraw()
histo <- histo +
  theme_linedraw(base_size = 8) +
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.margin = margin(2, 2, 2, 2),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
# Combined plot for EXP I Substrate effect on Growth
yr1_den_ht_plot_inset <- yr1_den_ht_plot +
  inset_element(
    histo,
    left = 0.7,    # x-position of left edge (0–1)
    bottom = 0.7,  # y-position of bottom edge (0–1)
    right = 0.98,   # x-position of right edge
    top = 0.98,     # y-position of top edge
    align_to = "panel")


# YEAR 2 - ROUNDNESS
# yr2.ht <- lmer(Height ~ Substrate + Flow_m + avg_Temp_c + Length + pre_Height + pre_Length + (1|Basket), data=yr2_df, REML = T)
# saveRDS(yr2.ht,file='data/yr2.ht.data')
check_model(yr2.ht, check = c("linearity", "homogeneity", "qq", "normality"))
summary(yr2.ht)
vc <- as.data.frame(VarCorr(yr2.ht))
var_basket <- vc$vcov[vc$grp == "Basket"]
var_resid <- vc$vcov[vc$grp == "Residual"]
var_basket / (var_basket + var_resid)
# yr2.ht_boot_fixef_lmer <- bootMer(yr2.ht, fixef_fun,nsim = 10000, use.u = FALSE, type = "parametric")
# saveRDS(yr2.ht_boot_fixef_lmer,file='data/yr2.ht_boot_fixef_lmer.data')
yr2.ht_ci_lmer <- apply(
  yr2.ht_boot_fixef_lmer$t,
  2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
yr2.ht_ci_lmer <- t(yr2.ht_ci_lmer)
colnames(yr2.ht_ci_lmer) <- c("CI_lower", "CI_upper")
data.frame(
  Estimate = fixef(yr2.ht),
  CI_lower = yr2.ht_ci_lmer[, "CI_lower"],
  CI_upper = yr2.ht_ci_lmer[, "CI_upper"])
# Flow and temperature have creditable effects!


# # YEAR 2 - ROUNDNESS - FLOW PLOTTING
# New data - keeping other fixed effects constant
flow_range <- seq(
  min(yr2_df$Flow_m, na.rm = TRUE),
  max(yr2_df$Flow_m, na.rm = TRUE),
  length.out = 100)
pred_flow_ht_data <- expand.grid(
  Flow_m = flow_range,
  Substrate = levels(factor(yr2_df$Substrate))[1],
  Length = mean(yr2_df$Length, na.rm = TRUE),
  pre_Length = mean(yr2_df$pre_Length, na.rm = TRUE),
  pre_Height = mean(yr2_df$pre_Height, na.rm = TRUE),
  avg_Temp_c = mean(yr2_df$avg_Temp_c, na.rm = TRUE))
# Predictions
pred_flow_ht_data$fit <- predict(yr2.ht, newdata = pred_flow_ht_data, re.form = NA)
# Bootstrap predictions
pred_flow_ht_FUN <- function(fit) {
  predict(fit, newdata = pred_flow_ht_data, re.form = NA)
}
# Run parametric bootstrapping
# pred_flow_ht_BOOT <- bootMer(yr2.ht, pred_flow_ht_FUN,
#                                  nsim = 10000, use.u = FALSE, type = "parametric")
# saveRDS(pred_flow_ht_BOOT,file='data/pred_flow_ht_BOOT.data')
# Lower bound of 95% CI
pred_flow_ht_data$lower <- apply(pred_flow_ht_BOOT$t, 2,
                                     quantile, probs = 0.025, na.rm = TRUE)
# Upper bound of 95% CI
pred_flow_ht_data$upper <- apply(pred_flow_ht_BOOT$t, 2,
                                     quantile, probs = 0.975, na.rm = TRUE)
boot_df <- data.frame(value = as.vector(pred_flow_ht_BOOT$t),
                      Flow_m = rep(pred_flow_ht_data$Flow_m,
                                   each = nrow(pred_flow_ht_BOOT$t)),
                      bootstrap_sample = paste0('bootstrap_',
                                                row_number(pred_flow_ht_BOOT$t)))
yr1_flow_ht_plot <- ggplot() +
  geom_line(data = boot_df, 
            aes(x = Flow_m, y = value, 
                group = bootstrap_sample),
            alpha = 0.05, size = 0.3, color='royalblue') +
  geom_line(data = pred_flow_ht_data,
            aes(x = Flow_m, y = fit),
            color = 'black',
            size = 1.5) +
  scale_x_continuous(limits = range(boot_df$Flow_m),
                     expand = c(0, 0)) +
  ylim(c(43.65,48.5)) +
  labs(x = "Flow velocity (m/s)",
       y = "Predicted shell height (mm)") +
  theme_linedraw() +
  theme(plot.margin = margin(2, 2, 2, 2))
max(pred_flow_ht_data$fit)-min(pred_flow_ht_data$fit)
# Get Histo for plot above
flow_boot <- yr2.ht_boot_fixef_lmer$t[, "Flow_m"]
flow_ci <- quantile(flow_boot, probs = c(0.025, 0.975), na.rm = TRUE)
histo <- ggplot(data.frame(flow_boot), aes(x = flow_boot)) +
  geom_density(fill = "grey80") +
  # xlim(c(-2,54)) +
  # ylim(c(0,.3)) +
  scale_y_continuous(limits=c(0, 1.32), breaks = seq(0, 1.3, by = .65)) +
  scale_x_continuous(limits=c(-15, 2), breaks = c(-14, -7, 0)) +
  geom_vline(xintercept = 0, linewidth = 0.7,
             linetype = "solid",color= "red") +
  geom_vline(xintercept = flow_ci,linewidth = 0.5,
             linetype = "solid",color= "blue",alpha=0.65) +
  labs(x = "Bootstrap estimates of Density effect",y = "Density") +
  theme_linedraw()
histo <- histo +
  theme_linedraw(base_size = 8) +
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.margin = margin(2, 2, 2, 2),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
# Combined plot for EXP II FLOW effect on ROUNDNESS
yr1_flow_ht_plot_inset <- yr1_flow_ht_plot +
  inset_element(
    histo,
    left = 0.7,    # x-position of left edge (0–1)
    bottom = 0.7,  # y-position of bottom edge (0–1)
    right = 0.98,   # x-position of right edge
    top = 0.98,     # y-position of top edge
    align_to = "panel")


# # YEAR 2 - ROUNDNESS - TEMP PLOTTING
# New data - keeping other fixed effects constant
temp_range <- seq(
  min(yr2_df$avg_Temp_c, na.rm = TRUE),
  max(yr2_df$avg_Temp_c, na.rm = TRUE),
  length.out = 100)
pred_temp_ht_data <- expand.grid(
  avg_Temp_c = temp_range,
  Substrate = levels(factor(yr2_df$Substrate))[1],
  Length = mean(yr2_df$Length, na.rm = TRUE),
  pre_Length = mean(yr2_df$pre_Length, na.rm = TRUE),
  pre_Height = mean(yr2_df$pre_Height, na.rm = TRUE),
  Flow_m = mean(yr2_df$Flow_m, na.rm = TRUE))
# Predictions
pred_temp_ht_data$fit <- predict(yr2.ht, newdata = pred_temp_ht_data, re.form = NA)
# Bootstrap predictions
pred_temp_ht_FUN <- function(fit) {
  predict(fit, newdata = pred_temp_ht_data, re.form = NA)
}
# Run parametric bootstrapping
# pred_temp_ht_BOOT <- bootMer(yr2.ht, pred_temp_ht_FUN,
#                                  nsim = 10000, use.u = FALSE, type = "parametric")
# saveRDS(pred_temp_ht_BOOT,file='data/pred_temp_ht_BOOT.data')
# Lower bound of 95% CI
pred_temp_ht_data$lower <- apply(pred_temp_ht_BOOT$t, 2,
                                     quantile, probs = 0.025, na.rm = TRUE)
# Upper bound of 95% CI
pred_temp_ht_data$upper <- apply(pred_temp_ht_BOOT$t, 2,
                                     quantile, probs = 0.975, na.rm = TRUE)
boot_df <- data.frame(value = as.vector(pred_temp_ht_BOOT$t),
                      avg_Temp_c = rep(pred_temp_ht_data$avg_Temp_c,
                                       each = nrow(pred_temp_ht_BOOT$t)),
                      bootstrap_sample = paste0('bootstrap_',
                                                row_number(pred_temp_ht_BOOT$t)))
yr1_temp_ht_plot <- ggplot() +
  geom_line(data = boot_df, 
            aes(x = avg_Temp_c, y = value, 
                group = bootstrap_sample),
            alpha = 0.05, size = 0.3, color='#DC143C') +
  geom_line(data = pred_temp_ht_data,
            aes(x = avg_Temp_c, y = fit),
            color = 'black',
            size = 1.5) +
  scale_x_continuous(limits = range(boot_df$avg_Temp_c),
                     expand = c(0, 0)) +
  ylim(c(43.65,48.5)) +
  labs(x = "Water temperature (°C)",
       y = "Predicted shell height (mm)") +
  theme_linedraw() +
  theme(plot.margin = margin(2, 2, 2, 2))
max(pred_temp_ht_data$fit)-min(pred_temp_ht_data$fit)
# Get Histo for plot above
temp_boot <- yr2.ht_boot_fixef_lmer$t[, "avg_Temp_c"]
temp_ci <- quantile(temp_boot, probs = c(0.025, 0.975), na.rm = TRUE)
histo <- ggplot(data.frame(temp_boot), aes(x = temp_boot)) +
  geom_density(fill = "grey80") +
  # xlim(c(-2,54)) +
  # ylim(c(0,.3)) +
  scale_y_continuous(limits=c(0, 1.32), breaks = seq(0, 1.3, by = .65)) +
  scale_x_continuous(limits=c(-15, 2), breaks = c(-14, -7, 0)) +
  geom_vline(xintercept = 0, linewidth = 0.7,
             linetype = "solid",color= "red") +
  geom_vline(xintercept = temp_ci,linewidth = 0.5,
             linetype = "solid",color= "blue",alpha=0.65) +
  labs(x = "Bootstrap estimates of Temp effect",y = "Density") +
  theme_linedraw()
histo <- histo +
  theme_linedraw(base_size = 8) +
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.margin = margin(2, 2, 2, 2),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
# Combined plot for EXP II TEMP effect on Growth
yr1_temp_ht_plot_inset <- yr1_temp_ht_plot +
  inset_element(
    histo,
    left = 0.7,    # x-position of left edge (0–1)
    bottom = 0.7,  # y-position of bottom edge (0–1)
    right = 0.98,   # x-position of right edge
    top = 0.98,     # y-position of top edge
    align_to = "panel")


left_panel <- (yr1_den_ht_plot_inset / plot_spacer() / plot_spacer()
               ) + 
  plot_layout(heights = c(1,0.05,1))
right_panel <- (yr1_flow_ht_plot_inset / plot_spacer() /
                  yr1_temp_ht_plot_inset) + 
  plot_layout(heights = c(1,0.05,1))
exp1_title <- textGrob(
  "Experiment I",
  gp = gpar(family = 'serif',fontface = "italic", fontsize = 11))
exp2_title <- textGrob(
  "Experiment II",
  gp = gpar(family = 'serif',fontface = "italic", fontsize = 11))
roundness_figure <- ((wrap_elements(exp1_title) |
                     wrap_elements(exp2_title)) /
                    ((left_panel | right_panel))) +
  plot_layout(
    heights = c(0.05, 1), widths = c(1,0.1,1)) 

roundness_figure <- roundness_figure &
  theme(text = element_text(family = "serif"))

ggsave(roundness_figure,filename='roundness_figure.png',path = paste0(path_figs),
       width = 6, height = 5.5, device='png', dpi=700)



###################################
# END: SHELL HEIGHT - ROUNDNESS
###################################




###################################
# SHELL WIDTH - INFLATION
###################################

# YEAR 1 - INFLATION
# yr1.wdt <- lmer(width ~ substrate + dist_m  + density + length + pre_width + pre_length + (1|basket), data=yr1_df, REML = T)
# saveRDS(yr1.wdt,file='data/yr1.wdt.data')
check_model(yr1.wdt, check = c("linearity", "homogeneity", "qq", "normality"))
summary(yr1.wdt)
vc <- as.data.frame(VarCorr(yr1.wdt))
var_basket <- vc$vcov[vc$grp == "basket"]
var_resid <- vc$vcov[vc$grp == "Residual"]
var_basket / (var_basket + var_resid)
# yr1.wdt_boot_fixef_lmer <- bootMer(yr1.wdt, fixef_fun,nsim = 10000, use.u = FALSE, type = "parametric")
# saveRDS(yr1.wdt_boot_fixef_lmer,file='data/yr1.wdt_boot_fixef_lmer.data')
yr1.wdt_ci_lmer <- apply(
  yr1.wdt_boot_fixef_lmer$t,
  2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
yr1.wdt_ci_lmer <- t(yr1.wdt_ci_lmer)
colnames(yr1.wdt_ci_lmer) <- c("CI_lower", "CI_upper")
data.frame(
  Estimate = fixef(yr1.wdt),
  CI_lower = yr1.wdt_ci_lmer[, "CI_lower"],
  CI_upper = yr1.wdt_ci_lmer[, "CI_upper"])
#Substrate & Density creditable effects!!



# # YEAR 1 - INFLATION - SUBSTRATE PLOTTING
# New data - keeping other fixed effects constant
pred_sub_wt_data <- expand.grid(
  substrate = levels(factor(yr1_df$substrate)),
  dist_m = mean(yr1_df$dist_m, na.rm = TRUE),
  density = levels(factor(yr1_df$density))[1],
  length = mean(yr1_df$length, na.rm = TRUE),
  pre_width = mean(yr1_df$pre_width, na.rm = TRUE),
  pre_length = mean(yr1_df$pre_length, na.rm = TRUE))
# Predictions
pred_sub_wt_data$fit <- predict(yr1.wdt, newdata = pred_sub_wt_data, re.form = NA)
# Bootstrap predictions
pred_sub_wt_FUN <- function(fit) {
  predict(fit, newdata = pred_sub_wt_data, re.form = NA)
}
# Run parametric bootstrapping
# pred_sub_wt_BOOT <- bootMer(yr1.wdt, pred_sub_wt_FUN,
#                                 nsim = 10000, use.u = FALSE, type = "parametric")
# saveRDS(pred_sub_wt_BOOT,file='data/pred_sub_wt_BOOT.data')
# Lower bound of 95% CI
pred_sub_wt_data$lower <- apply(pred_sub_wt_BOOT$t, 2,
                                    quantile, probs = 0.025, na.rm = TRUE)
# Upper bound of 95% CI
pred_sub_wt_data$upper <- apply(pred_sub_wt_BOOT$t, 2,
                                    quantile, probs = 0.975, na.rm = TRUE)
boot_df <- data.frame(value = as.vector(pred_sub_wt_BOOT$t),
                      substrate = rep(pred_sub_wt_data$substrate, each = nrow(pred_sub_wt_BOOT$t)))
yr1_sub_wt_plot <- ggplot() +
  geom_violin(data = boot_df, 
              aes(x = substrate, y = value, fill = substrate), trim = FALSE) +
  scale_fill_manual(values = c('#AF804F','#7F5112')) +
  geom_point(data = pred_sub_wt_data,
             aes(x = substrate, y = fit), size = 3) +
  geom_errorbar(data = pred_sub_wt_data,
                aes(x = substrate, ymin = lower, ymax = upper), width = 0.15) +
  scale_x_discrete(labels= c('Sand', 'Gravel')) +
  ylim(c(23.3,25.25)) +
  labs(x = "Substrate size", y = "Predicted shell width (mm)") + 
  guides(fill="none") +
  theme_linedraw() +
  theme(plot.margin = margin(2, 2, 2, 2))
# Get Histo for plot above
sub_boot <- yr1.wdt_boot_fixef_lmer$t[, "substrategravel"]
sub_ci <- quantile(sub_boot, probs = c(0.025, 0.975), na.rm = TRUE)
histo <- ggplot(data.frame(sub_boot), aes(x = sub_boot)) +
  geom_density(fill = "grey80") +
  # xlim(c(-20,2)) +
  # ylim(c(0,.2)) +
  scale_y_continuous(limits=c(0, 2.8), breaks = seq(0, 2.8, by = 1.4)) +
  scale_x_continuous(limits=c(-1.3, 0.25), breaks = c(-1.2, -.6, 0)) +
  geom_vline(xintercept = 0, linewidth = 0.7,
             linetype = "solid",color= "red") +
  geom_vline(xintercept = sub_ci,linewidth = 0.5,
             linetype = "solid",color= "blue",alpha=0.65) +
  labs(x = "Bootstrap estimates of Substrate effect",y = "Density") +
  theme_linedraw()
histo <- histo +
  theme_linedraw(base_size = 8) +
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.margin = margin(2, 2, 2, 2),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
# Combined plot for EXP I Substrate effect on Growth
yr1_sub_wt_plot_inset <- yr1_sub_wt_plot +
  inset_element(
    histo,
    left = 0.7,    # x-position of left edge (0–1)
    bottom = 0.7,  # y-position of bottom edge (0–1)
    right = 0.98,   # x-position of right edge
    top = 0.98,     # y-position of top edge
    align_to = "panel")



# # YEAR 1 - INFLATION - DENSITY PLOTTING
# New data - keeping other fixed effects constant
pred_den_wt_data <- expand.grid(
  substrate = levels(factor(yr1_df$substrate))[1],
  dist_m = mean(yr1_df$dist_m, na.rm = TRUE),
  density = levels(factor(yr1_df$density)),
  length = mean(yr1_df$length, na.rm = TRUE),
  pre_width = mean(yr1_df$pre_width, na.rm = TRUE),
  pre_length = mean(yr1_df$pre_length, na.rm = TRUE))
# Predictions
pred_den_wt_data$fit <- predict(yr1.wdt, newdata = pred_den_wt_data, re.form = NA)
# Bootstrap predictions
pred_den_wt_FUN <- function(fit) {
  predict(fit, newdata = pred_den_wt_data, re.form = NA)
}
# Run parametric bootstrapping
# pred_den_wt_BOOT <- bootMer(yr1.wdt, pred_den_wt_FUN,
#                                 nsim = 10000, use.u = FALSE, type = "parametric")
# saveRDS(pred_den_wt_BOOT,file='data/pred_den_wt_BOOT.data')
# Lower bound of 95% CI
pred_den_wt_data$lower <- apply(pred_den_wt_BOOT$t, 2,
                                    quantile, probs = 0.025, na.rm = TRUE)
# Upper bound of 95% CI
pred_den_wt_data$upper <- apply(pred_den_wt_BOOT$t, 2,
                                    quantile, probs = 0.975, na.rm = TRUE)
boot_df <- data.frame(value = as.vector(pred_den_wt_BOOT$t),
                      density = rep(pred_den_wt_data$density,
                                    each = nrow(pred_den_wt_BOOT$t)))
yr1_den_wt_plot <- ggplot() +
  geom_violin(data = boot_df, 
              aes(x = density, y = value, fill = density), trim = FALSE) +
  scale_fill_manual(values = c('gray90','gray40')) +
  geom_point(data = pred_den_wt_data,
             aes(x = density, y = fit), size = 3) +
  geom_errorbar(data = pred_den_wt_data,
                aes(x = density, ymin = lower, ymax = upper), width = 0.15) +
  scale_x_discrete(labels= c('Low - 10 indivs.', 'High - 50 indivs.')) +
  ylim(c(23.3,25.25)) +
  labs(x = "Mussel density", y = "Predicted shell width (mm)") + 
  guides(fill="none") +
  theme_linedraw() +
  theme(plot.margin = margin(2, 2, 2, 2))
# Get Histo for plot above
den_boot <- yr1.wdt_boot_fixef_lmer$t[, "densityhigh"]
den_ci <- quantile(den_boot, probs = c(0.025, 0.975), na.rm = TRUE)
histo <- ggplot(data.frame(den_boot), aes(x = den_boot)) +
  geom_density(fill = "grey80") +
  # xlim(c(-20,2)) +
  # ylim(c(0,.2)) +
  scale_y_continuous(limits=c(0, 2.8), breaks = seq(0, 2.8, by = 1.4)) +
  scale_x_continuous(limits=c(-1.3, 0.25), breaks = c(-1.2, -.6, 0)) +
  geom_vline(xintercept = 0, linewidth = 0.7,
             linetype = "solid",color= "red") +
  geom_vline(xintercept = den_ci,linewidth = 0.5,
             linetype = "solid",color= "blue",alpha=0.65) +
  labs(x = "Bootstrap estimates of Density effect",y = "Density") +
  theme_linedraw()
histo <- histo +
  theme_linedraw(base_size = 8) +
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.margin = margin(2, 2, 2, 2),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
# Combined plot for EXP I Substrate effect on Growth
yr1_den_wt_plot_inset <- yr1_den_wt_plot +
  inset_element(
    histo,
    left = 0.7,    # x-position of left edge (0–1)
    bottom = 0.7,  # y-position of bottom edge (0–1)
    right = 0.98,   # x-position of right edge
    top = 0.98,     # y-position of top edge
    align_to = "panel")



# YEAR 2 - INFLATION
# yr2.wdt <- lmer(Width ~ Substrate + Flow_m + avg_Temp_c + Length + pre_Width + pre_Length + (1|Basket), data=yr2_df, REML = T)
# saveRDS(yr2.wdt,file='data/yr2.wdt.data')
check_model(yr2.wdt, check = c("linearity", "homogeneity", "qq", "normality"))
summary(yr2.wdt)
vc <- as.data.frame(VarCorr(yr2.wdt))
var_basket <- vc$vcov[vc$grp == "Basket"]
var_resid <- vc$vcov[vc$grp == "Residual"]
var_basket / (var_basket + var_resid)
# yr2.wdt_boot_fixef_lmer <- bootMer(yr2.wdt, fixef_fun,nsim = 10000, use.u = FALSE, type = "parametric")
# saveRDS(yr2.wdt_boot_fixef_lmer,file='data/yr2.wdt_boot_fixef_lmer.data')
yr2.wdt_ci_lmer <- apply(
  yr2.wdt_boot_fixef_lmer$t,
  2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
yr2.wdt_ci_lmer <- t(yr2.wdt_ci_lmer)
colnames(yr2.wdt_ci_lmer) <- c("CI_lower", "CI_upper")
data.frame(
  Estimate = fixef(yr2.wdt),
  CI_lower = yr2.wdt_ci_lmer[, "CI_lower"],
  CI_upper = yr2.wdt_ci_lmer[, "CI_upper"])
# No creditable effects!!


left_panel <- (yr1_sub_wt_plot_inset / plot_spacer() / yr1_den_wt_plot_inset
) + 
  plot_layout(heights = c(1,0.05,1))
exp1_title <- textGrob(
  "Experiment I",
  gp = gpar(family = 'serif',fontface = "italic", fontsize = 11))
inflation_figure <- (wrap_elements(exp1_title)  /
                       (left_panel )) +
  plot_layout(
    heights = c(0.05, 1), widths = c(1)) 

inflation_figure <- inflation_figure &
  theme(text = element_text(family = "serif"))

ggsave(inflation_figure,filename='inflation_figure.png',path = paste0(path_figs),
       width = 6/2, height = 5.5, device='png', dpi=700)



###################################
# END: SHELL WIDTH - INFLATION
###################################




###################################
# SHELL THICKNESS INDEX
###################################

# YEAR 2 - SHELL THICKNESS - MODEL
# yr2.sti <- lmer(STI ~ Substrate + Flow_m + avg_Temp_c + Length + (1|Basket), data=yr2_df, REML = T)
# saveRDS(yr2.sti,file='data/yr2.sti.data')
check_model(yr2.sti, check = c("linearity", "homogeneity", "qq", "normality"))
summary(yr2.sti)
vc <- as.data.frame(VarCorr(yr2.sti))
var_basket <- vc$vcov[vc$grp == "Basket"]
var_resid <- vc$vcov[vc$grp == "Residual"]
var_basket / (var_basket + var_resid)
# yr2.sti_boot_fixef_lmer <- bootMer(yr2.sti, fixef_fun,nsim = 10000, use.u = FALSE, type = "parametric")
# saveRDS(yr2.sti_boot_fixef_lmer,file='data/yr2.sti_boot_fixef_lmer.data')
yr2.sti_ci_lmer <- apply(
  yr2.sti_boot_fixef_lmer$t,
  2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
yr2.sti_ci_lmer <- t(yr2.sti_ci_lmer)
colnames(yr2.sti_ci_lmer) <- c("CI_lower", "CI_upper")
data.frame(
  Estimate = fixef(yr2.sti),
  CI_lower = yr2.sti_ci_lmer[, "CI_lower"],
  CI_upper = yr2.sti_ci_lmer[, "CI_upper"])
# Flow and Temp have creditable effects!!



# # YEAR 2 - SHELL THICKNESS INDEX - FLOW PLOTTING
# New data - keeping other fixed effects constant
flow_range <- seq(
  min(yr2_df$Flow_m, na.rm = TRUE),
  max(yr2_df$Flow_m, na.rm = TRUE),
  length.out = 100)
pred_flow_sti_data <- expand.grid(
  Flow_m = flow_range,
  Substrate = levels(factor(yr2_df$Substrate))[1],
  avg_Temp_c = mean(yr2_df$avg_Temp_c, na.rm = TRUE),
  Length = mean(yr2_df$Length, na.rm = TRUE))
# Predictions
pred_flow_sti_data$fit <- predict(yr2.sti, newdata = pred_flow_sti_data, re.form = NA)
# Bootstrap predictions
pred_flow_sti_FUN <- function(fit) {
  predict(fit, newdata = pred_flow_sti_data, re.form = NA)
}
# Run parametric bootstrapping
# pred_flow_sti_BOOT <- bootMer(yr2.sti, pred_flow_sti_FUN,
#                              nsim = 10000, use.u = FALSE, type = "parametric")
# saveRDS(pred_flow_sti_BOOT,file='data/pred_flow_sti_BOOT.data')
# Lower bound of 95% CI
pred_flow_sti_data$lower <- apply(pred_flow_sti_BOOT$t, 2,
                                 quantile, probs = 0.025, na.rm = TRUE)
# Upper bound of 95% CI
pred_flow_sti_data$upper <- apply(pred_flow_sti_BOOT$t, 2,
                                 quantile, probs = 0.975, na.rm = TRUE)
boot_df <- data.frame(value = as.vector(pred_flow_sti_BOOT$t),
                      Flow_m = rep(pred_flow_sti_data$Flow_m,
                                   each = nrow(pred_flow_sti_BOOT$t)),
                      bootstrap_sample = paste0('bootstrap_',
                                                row_number(pred_flow_sti_BOOT$t)))
yr1_flow_sti_plot <- ggplot() +
  geom_line(data = boot_df, 
            aes(x = Flow_m, y = value, 
                group = bootstrap_sample),
            alpha = 0.05, size = 0.3, color='royalblue') +
  geom_line(data = pred_flow_sti_data,
            aes(x = Flow_m, y = fit),
            color = 'black',
            size = 1.5) +
  scale_x_continuous(limits = range(boot_df$Flow_m),
                     expand = c(0, 0)) +
  ylim(c(0.9,1.63)) +
  labs(x = "Flow velocity (m/s)",
       y = "Predicted shell thickness index") +
  theme_linedraw() +
  theme(plot.margin = margin(2, 2, 2, 2))
# Get Histo for plot above
flow_boot <- yr2.sti_boot_fixef_lmer$t[, "Flow_m"]
flow_ci <- quantile(flow_boot, probs = c(0.025, 0.975), na.rm = TRUE)
histo <- ggplot(data.frame(flow_boot), aes(x = flow_boot)) +
  geom_density(fill = "grey80") +
  # xlim(c(-2,54)) +
  # ylim(c(0,.3)) +
  scale_y_continuous(limits=c(0, 12), breaks = seq(0, 12, by = 6)) +
  scale_x_continuous(limits=c(-0.6, 1.8), breaks = c(0, .8, 1.6)) +
  geom_vline(xintercept = 0, linewidth = 0.7,
             linetype = "solid",color= "red") +
  geom_vline(xintercept = flow_ci,linewidth = 0.5,
             linetype = "solid",color= "blue",alpha=0.65) +
  labs(x = "Bootstrap estimates of Density effect",y = "Density") +
  theme_linedraw()
histo <- histo +
  theme_linedraw(base_size = 8) +
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.margin = margin(2, 2, 2, 2),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
# Combined plot for EXP II FLOW effect on ROUNDNESS
yr1_flow_sti_plot_inset <- yr1_flow_sti_plot +
  inset_element(
    histo,
    left = 0.02,    # x-position of left edge (0–1)
    bottom = 0.7,  # y-position of bottom edge (0–1)
    right = 0.3,   # x-position of right edge
    top = 0.98,     # y-position of top edge
    align_to = "panel")


# # YEAR 2 - SHELL THICKNESS INDEX - TEMP PLOTTING
# New data - keeping other fixed effects constant
temp_range <- seq(
  min(yr2_df$avg_Temp_c, na.rm = TRUE),
  max(yr2_df$avg_Temp_c, na.rm = TRUE),
  length.out = 100)
pred_temp_sti_data <- expand.grid(
  avg_Temp_c = temp_range,
  Substrate = levels(factor(yr2_df$Substrate))[1],
  Flow_m = mean(yr2_df$Flow_m, na.rm = TRUE),
  Length = mean(yr2_df$Length, na.rm = TRUE))
# Predictions
pred_temp_sti_data$fit <- predict(yr2.sti, newdata = pred_temp_sti_data, re.form = NA)
# Bootstrap predictions
pred_temp_sti_FUN <- function(fit) {
  predict(fit, newdata = pred_temp_sti_data, re.form = NA)
}
# Run parametric bootstrapping
# pred_temp_sti_BOOT <- bootMer(yr2.sti, pred_temp_sti_FUN,
#                              nsim = 10000, use.u = FALSE, type = "parametric")
# saveRDS(pred_temp_sti_BOOT,file='data/pred_temp_sti_BOOT.data')
# Lower bound of 95% CI
pred_temp_sti_data$lower <- apply(pred_temp_sti_BOOT$t, 2,
                                 quantile, probs = 0.025, na.rm = TRUE)
# Upper bound of 95% CI
pred_temp_sti_data$upper <- apply(pred_temp_sti_BOOT$t, 2,
                                 quantile, probs = 0.975, na.rm = TRUE)
boot_df <- data.frame(value = as.vector(pred_temp_sti_BOOT$t),
                      avg_Temp_c = rep(pred_temp_sti_data$avg_Temp_c,
                                       each = nrow(pred_temp_sti_BOOT$t)),
                      bootstrap_sample = paste0('bootstrap_',
                                                row_number(pred_temp_sti_BOOT$t)))
yr1_temp_sti_plot <- ggplot() +
  geom_line(data = boot_df, 
            aes(x = avg_Temp_c, y = value, 
                group = bootstrap_sample),
            alpha = 0.05, size = 0.3, color='#DC143C') +
  geom_line(data = pred_temp_sti_data,
            aes(x = avg_Temp_c, y = fit),
            color = 'black',
            size = 1.5) +
  scale_x_continuous(limits = range(boot_df$avg_Temp_c),
                     expand = c(0, 0)) +
  ylim(c(0.9,1.63)) +
  labs(x = "Water temperature (°C)",
       y = "Predicted shell thickness index") +
  theme_linedraw() +
  theme(plot.margin = margin(2, 2, 2, 2))
# Get Histo for plot above
temp_boot <- yr2.sti_boot_fixef_lmer$t[, "avg_Temp_c"]
temp_ci <- quantile(temp_boot, probs = c(0.025, 0.975), na.rm = TRUE)
histo <- ggplot(data.frame(temp_boot), aes(x = temp_boot)) +
  geom_density(fill = "grey80") +
  # xlim(c(-2,54)) +
  # ylim(c(0,.3)) +
  scale_y_continuous(limits=c(0, 12), breaks = seq(0, 12, by = 6)) +
  scale_x_continuous(limits=c(-0.6, 1.8), breaks = c(0, .8, 1.6)) +
  geom_vline(xintercept = 0, linewidth = 0.7,
             linetype = "solid",color= "red") +
  geom_vline(xintercept = temp_ci,linewidth = 0.5,
             linetype = "solid",color= "blue",alpha=0.65) +
  labs(x = "Bootstrap estimates of Temp effect",y = "Density") +
  theme_linedraw()
histo <- histo +
  theme_linedraw(base_size = 8) +
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.margin = margin(2, 2, 2, 2),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
# Combined plot for EXP II TEMP effect on Growth
yr1_temp_sti_plot_inset <- yr1_temp_sti_plot +
  inset_element(
    histo,
    left = 0.02,    # x-position of left edge (0–1)
    bottom = 0.7,  # y-position of bottom edge (0–1)
    right = 0.3,   # x-position of right edge
    top = 0.98,     # y-position of top edge
    align_to = "panel")


right_panel <- (yr1_flow_sti_plot_inset / plot_spacer() / yr1_temp_sti_plot_inset
) + 
  plot_layout(heights = c(1,0.05,1))
exp2_title <- textGrob(
  "Experiment II",
  gp = gpar(family = 'serif',fontface = "italic", fontsize = 11))
sti_figure <- (wrap_elements(exp2_title)  /
                       (right_panel )) +
  plot_layout(
    heights = c(0.05, 1), widths = c(1)) 

sti_figure <- sti_figure &
  theme(text = element_text(family = "serif"))

ggsave(sti_figure,filename='sti_figure.png',path = paste0(path_figs),
       width = 6/2, height = 5.5, device='png', dpi=700)




###################################
# END: SHELL THICKNESS INDEX
###################################




###################################
###################################
# START: LINEAR REGRESSIONS & BOOTSTRAPPING - SENSITIVITY ANALYSIS
###################################
###################################

###################################
# SHELL LENGTH - GROWTH RATE
###################################

# YEAR 1 - GROWTH
yr1.len <- lm(
  length_growth ~ substrate + dist_m + density + pre_length, 
  data=yr1_df)
check_model(yr1.len, check = c("linearity", "homogeneity", "qq", "normality"))
lm_fixef_fun <- function(data, indices) {
  fit <- lm(
    length_growth ~ substrate + dist_m + density + pre_length, 
    data = data[indices, ])
  coef(fit)}
boot_fixef_lm <- boot(data = yr1_df, statistic = lm_fixef_fun, R = 10000)
ci_lm <- apply(boot_fixef_lm$t,2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
ci_lm <- t(ci_lm)
colnames(ci_lm) <- c("CI_lower", "CI_upper")
data.frame(
  Estimate = coef(yr1.len),
  CI_lower = ci_lm[, "CI_lower"],
  CI_upper = ci_lm[, "CI_upper"])

# YEAR 2 - GROWTH
yr2.len <- lm(
  Length_growth ~ Substrate + Flow_m + avg_Temp_c + pre_Length,
  data=yr2_df)
check_model(yr2.len, check = c("linearity", "homogeneity", "qq", "normality"))
lm_fixef_fun <- function(data, indices) {
  fit <- lm(
    Length_growth ~ Substrate + Flow_m + avg_Temp_c + pre_Length,
    data = data[indices, ])
  coef(fit)}
boot_fixef_lm <- boot(data = yr2_df, statistic = lm_fixef_fun, R = 10000)
ci_lm <- apply(boot_fixef_lm$t,2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
ci_lm <- t(ci_lm)
colnames(ci_lm) <- c("CI_lower", "CI_upper")
data.frame(
  Estimate = coef(yr2.len),
  CI_lower = ci_lm[, "CI_lower"],
  CI_upper = ci_lm[, "CI_upper"])


###################################
# SHELL HEIGHT - ROUNDNESS
###################################

# YEAR 1 - ROUNDNESS
yr1.ht <- lm(
  height ~ substrate + dist_m  + density + length + pre_height + pre_length,
  data=yr1_df)
check_model(yr1.ht, check = c("linearity", "homogeneity", "qq", "normality"))
lm_fixef_fun <- function(data, indices) {
  fit <- lm(
    height ~ substrate + dist_m  + density + length + pre_height + pre_length,
    data = data[indices, ])
  coef(fit)}
boot_fixef_lm <- boot(data = yr1_df, statistic = lm_fixef_fun, R = 10000)
ci_lm <- apply(boot_fixef_lm$t,2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
ci_lm <- t(ci_lm)
colnames(ci_lm) <- c("CI_lower", "CI_upper")
data.frame(
  Estimate = coef(yr1.ht),
  CI_lower = ci_lm[, "CI_lower"],
  CI_upper = ci_lm[, "CI_upper"])

# YEAR 2 - ROUNDNESS
yr2.ht <- lm(
  Height ~ Substrate + Flow_m + avg_Temp_c + Length + pre_Height + pre_Length,
  data=yr2_df)
check_model(yr2.ht, check = c("linearity", "homogeneity", "qq", "normality"))
lm_fixef_fun <- function(data, indices) {
  fit <- lm(
    Height ~ Substrate + Flow_m + avg_Temp_c + Length + pre_Height + pre_Length,
    data = data[indices, ])
  coef(fit)}
boot_fixef_lm <- boot(data = yr2_df, statistic = lm_fixef_fun, R = 10000)
ci_lm <- apply(boot_fixef_lm$t,2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
ci_lm <- t(ci_lm)
colnames(ci_lm) <- c("CI_lower", "CI_upper")
data.frame(
  Estimate = coef(yr2.ht),
  CI_lower = ci_lm[, "CI_lower"],
  CI_upper = ci_lm[, "CI_upper"])


# YEAR 1 - INFLATION
yr1.wdt <- lm(
  width ~ substrate + dist_m  + density + length + pre_width + pre_length,
  data=yr1_df)
check_model(yr1.wdt, check = c("linearity", "homogeneity", "qq", "normality"))
lm_fixef_fun <- function(data, indices) {
  fit <- lm(
    width ~ substrate + dist_m  + density + length + pre_width + pre_length,
    data = data[indices, ])
  coef(fit)}
boot_fixef_lm <- boot(data = yr1_df, statistic = lm_fixef_fun, R = 10000)
ci_lm <- apply(boot_fixef_lm$t,2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
ci_lm <- t(ci_lm)
colnames(ci_lm) <- c("CI_lower", "CI_upper")
data.frame(
  Estimate = coef(yr1.wdt),
  CI_lower = ci_lm[, "CI_lower"],
  CI_upper = ci_lm[, "CI_upper"])


# YEAR 2 - INFLATION
yr2.wdt <- lm(
  Width ~ Substrate + Flow_m + avg_Temp_c + Length + pre_Width + pre_Length,
  data=yr2_df)
check_model(yr2.wdt, check = c("linearity", "homogeneity", "qq", "normality"))
lm_fixef_fun <- function(data, indices) {
  fit <- lm(
    Width ~ Substrate + Flow_m + avg_Temp_c + Length + pre_Width + pre_Length,
    data = data[indices, ])
  coef(fit)}
boot_fixef_lm <- boot(data = yr2_df, statistic = lm_fixef_fun, R = 10000)
ci_lm <- apply(boot_fixef_lm$t,2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
ci_lm <- t(ci_lm)
colnames(ci_lm) <- c("CI_lower", "CI_upper")
data.frame(
  Estimate = coef(yr2.wdt),
  CI_lower = ci_lm[, "CI_lower"],
  CI_upper = ci_lm[, "CI_upper"])


###################################
# SHELL THICKNESS INDEX
###################################

# YEAR 2 - SHELL THICKNESS
yr2.sti <- lm(
  STI ~ Substrate + Flow_m + avg_Temp_c + Length,
  data=yr2_df)
check_model(yr2.sti, check = c("linearity", "homogeneity", "qq", "normality"))
lm_fixef_fun <- function(data, indices) {
  fit <- lm(
    STI ~ Substrate + Flow_m + avg_Temp_c + Length,
    data = data[indices, ])
  coef(fit)}
boot_fixef_lm <- boot(data = yr2_df, statistic = lm_fixef_fun, R = 10000)
ci_lm <- apply(boot_fixef_lm$t,2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
ci_lm <- t(ci_lm)
colnames(ci_lm) <- c("CI_lower", "CI_upper")
data.frame(
  Estimate = coef(yr2.sti),
  CI_lower = ci_lm[, "CI_lower"],
  CI_upper = ci_lm[, "CI_upper"])


###################################
# END: LINEAR REGRESSIONS & BOOTSTRAPPING - SENSITIVITY ANALYSIS
###################################


###################################
# START: CALCULATE BIOMASS OF DENSITY TREATMENTS
###################################


# Density settings reasonable?
mean(may$pre_Mass)
mean(may$pre_Length)
mean(yr1_df$pre_length)

df_new <- data.frame(mass = c(yr2_df$pre_Mass,yr2_df$Mass),
                     length= c(yr2_df$pre_Length,yr2_df$Length))

mod <- lm(mass ~ poly(length, 2), data = df_new)
summary(mod)
newx <- data.frame(length = seq(38, max(df_new$length) * 1.2, length.out = 200))
pred <- predict(mod, newx, se.fit = TRUE)
newx$y <- pred$fit

ggplot(df_new, aes(length, mass)) +
  geom_point() +
  geom_line(data = newx, aes(length, y)) +
  xlim(25,120)

# 41.7 mm
mean(yr1_df$pre_length)
# Pred shows that mussels at 41.63 mm are 4.61 g
newx

highden_grams <- 4.61*50
lowden_grams <- 4.61*10
meso_area <- .2*.2*3.14159
multiplier<- 1/meso_area
# 1834.3 g/m2
highden_grams*multiplier
# 366.9 g/m2
lowden_grams*multiplier



###################################
# END: KEOGH ET AL.
###################################
