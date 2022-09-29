# Set working directory

setwd("~/Desktop/Projects/Microbiome/Analysis/")

## Load all required packages

library(lme4)
library(lmerTest)



library(ggplot2)
library(phyloseq)
library(marmotdata)
library(MuMIn)
library(see)
library(performance)


## Read in the datasets
## First the full dataset with all ages included

all_ages <- read.csv("sample_metadata_allyears.csv")

family_all_ages <- read.csv("sample_metadata_allyears_family.csv")

bac_plot <- read.csv("sample_metadata_allyears_BAC.csv")

muri_plot <- read.csv("muri_plot_data.csv")

muri_plot_UV <- read.csv("muri_plot_data_UV.csv")

## Just juveniles

juv <- read.csv("juvenile_allyears.csv")

## Just yearlings

yr <- read.csv("yearlings_allyears.csv")

## Just adults

ad <- read.csv("adults_allyears.csv")

## Dataset with Bacteroidetes and up valley colonies

bac_uv <- read.csv("sample_metadata_bac_uv.csv")


## Check distribution of the response variable

hist(all_ages$mgr)
shapiro.test(all_ages$mgr) ## Normally distributed


## Fit linear mixed models
## First use all age dataset


## For the phylum Firmicutes 
## Valley position as fixed effect and year as a random effect
## Scale and center firmicutes

scale_Firmicutes <- scale(all_ages$Firmicutes, center = TRUE, scale = TRUE)

model_f1_all <- lmer(mgr ~ scale_Firmicutes + valley_position + (1|year) + (1|uid), data = all_ages)
summary(model_f1_all)
r.squaredGLMM(model_f1_all)
## Valley position as random effect and year as a random effect

model_f2_all <- lmer(mgr ~ Firmicutes + age_class + (1|valley_position) + (1|year) + (1|uid), data = all_ages)

## Valley position and year both as fixed effects

model_f3_all <- lmer(mgr ~ Firmicutes + valley_position + age_class + year + (1|uid), data = all_ages)

## Valley position as random effect and year as fixed effect

model_f4_all <- lmer(mgr ~ Firmicutes + age_class + year + (1|valley_position) + (1|uid), data = all_ages)

### Some interactions

model_I1_all <- lmer(mgr ~ valley_position + scale_Bacteroidetes*valley_position + (1|year) + (1|uid), data = all_ages)
summary(model_I1_all)
r.squaredGLMM(model_I1_all)
model_I2_all <- lmer(mgr ~ scale_Firmicutes + valley + age_c + (1|year) + (1|uid), data = all_ages)
hist(all_ages$mgr)
hist(all_ages$Firmicutes)

model_noF_all <- lmer(mgr ~ valley_position + age_class + (1|year) + (1|uid), data = all_ages)
## For the phylum of Bacteroidetes
## Valley position as fixed effect and year as a random effect
## Scale and center bacteroidetes

scale_Bacteroidetes <- scale(all_ages$Bacteroidetes, center = TRUE, scale = TRUE)

model_b1_all <- lmer(mgr ~ scale_Bacteroidetes + valley_position + age_class + (1|year) + (1|uid), data = all_ages)

## Valley position as random effect and year as a random effect

model_b2_all <- lmer(mgr ~ Bacteroidetes + age_class + (1|valley_position) + (1|year) + (1|uid), data = all_ages)

## Valley position and year both as fixed effects

model_b3_all <- lmer(mgr ~ Bacteroidetevalley_position + age_class + year + (1|uid), data = all_ages)

## Valley position as random effect and year as fixed effect

model_b4_all <- lmer(mgr ~ Bacteroidetes + age_class + year + (1|valley_position) + (1|uid), data = all_ages)



## Explore some interactions between age and bacterial groups 
## Functional class, age, and interaction between the two 

summary(model_f1_all)
summary(model_b1_all)
summary(model_I1_all)
anova(model_I1_all)
summary(model2_all)
summary(model3_all)
summary(model4_all)
summary(model5_all)

summary(model_c_all)
anova(model_c_all)

#### Models at the family level

scale_rumi <- scale(family_all_ages$rumi_fam, center = TRUE, scale = TRUE)

model_rumi_fam <- lmer(mgr ~ valley_position + age_class + (1|year) + (1|uid), data = all_ages)
summary(model_rumi_fam)


scale_muri <- scale(family_all_ages$muri_fam, center = TRUE, scale = TRUE)

model_muri_fam <- lmer(mgr ~scale_muri*valley_position + valley_position + age_class + (1|year) + (1|uid), data = all_ages)
summary(model_muri_fam)


## Check the residuals
check_model(model_rumi_fam)

## Check effect sizes

r.squaredGLMM(model_muri_fam)


##### Plotting

install.packages("sjPlot")
library(sjPlot)

## Graph and Table for Firmicutes model
## Scale and Center first 
scale_Firmicutes <- scale(all_ages$Firmicutes, center = TRUE, scale = TRUE)
ggplot(all_ages, aes(x = scale_Firmicutes, y = mgr)) + geom_smooth(method = "lm", se = FALSE) + geom_point()

tab_model(model_f1_all, 
          show.re.var= TRUE, 
          pred.labels =c("(Intercept)", "Firmicutes", "Valley Position", "Age [Juveniles]", "Age [Yearlings]"),
          dv.labels= "Effects of Firmicutes, Age, and Valley Position on Mass Gain Rate")


## Graph and Table for Bacteroidetes Model

## Scale and center first 

scale_Bacteroidetes <- scale(bac_plot$Bacteroidetes, center = TRUE, scale = TRUE)
ggplot(all_ages, aes(x = scale_Bacteroidetes, y = scale_mgr)) + geom_smooth(method = "lm", se = FALSE) + geom_point()

tab_model(model_b1_all, 
          show.re.var= TRUE, 
          pred.labels =c("(Intercept)", "Bacteroidetes", "Valley Position", "Age [Juveniles]", "Age [Yearlings]"),
          dv.labels= "Effects of Bacteroidetes, Age, and Valley Position on Mass Gain Rate")

## Plot of the interaction of Bacteroidetes and valley position

tab_model(model_I1_all, 
          show.re.var= TRUE, 
          pred.labels =c("(Intercept)", "Bacteroidetes", "Valley Position", "Age [Juveniles]", "Age [Yearlings]", "Bacteroidetes*[Up-Valley]"),
          dv.labels= "Effect of Bacteroidetes, Valley Position, and Age Class on Mass Gain Rate with Interaction")


## Plot of effect sizes

plot_model(all_phyla, 
           axis.labels=c("Firmicutes", "Bacteroidetes", "Valley Position", "Age [Juveniles]", "Age [Yearlings]"),
           show.values=TRUE, show.p=TRUE,
           title="Effect of Bacterial Phyla, Valley Position, and Age on Mass Gain Rate")

### Make some figures for the paper!
## First do Firmicutes and mass gain rate

library(ggplot2)
# Scatter plot representing the relationship between Firmicutes and Mass Gain Rate on the full dataset
firm_plot <- ggplot(all_ages, aes(x=scale_Firmicutes, y=mgr)) + geom_point() + geom_smooth(method=lm, se=TRUE) + xlab(expression(paste(italic("Firmicutes "), "Abundance"))) + ylab("Mass Gain Rate (g/day)")
firm_plot + theme_linedraw(base_size = 22)

## Plot for Bacteroidetes and up valley

bac_graph <- ggplot(bac_plot, aes(x = scale_Bacteroidetes, y = mgr, color = valley_position)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) + xlab(expression(paste(italic("Bacteroidetes "), "Abundance"))) + ylab("Mass Gain Rate (g/day)")
bac_graph + theme_linedraw(base_size = 22) + labs(color = "Valley Position")

  
  
install.packages("ggtext") 
library(ggtext)

library(marmotdata)


install.packages("ggpubr")
library(ggpubr)
install.packages("cowplot")
library(cowplot)


#### FINAL PLOTS ####

# Scatter plot representing the relationship between Firmicutes and Mass Gain Rate on the full dataset

firm_plot <- ggplot(all_ages, aes(x=scale_Firmicutes, y=mgr)) + geom_point() + geom_smooth(method = "lm", se = TRUE) +  xlab(expression(paste(italic("Firmicutes "), "Abundance"))) + ylab("Mass Gain Rate (g/day)")
firm_plot_final <- firm_plot + theme_linedraw(base_size = 22) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(breaks=c(-2,-1, 0, 1, 2, 3, 4))


# Scatter plot representing the relationship between Bacteroidetes and Mass Gain Rate on the full dataset

bac_graph <- ggplot(bac_plot, aes(x = scale_Bacteroidetes, y = mgr, color = valley_position)) +
  geom_point() +
  geom_smooth(aes(fill = valley_position), method = "lm", se = TRUE) + xlab(expression(paste(italic("Bacteroidetes "), "Abundance"))) + ylab("")
bac_plot_final <- bac_graph + theme_linedraw(base_size = 22) + theme(legend.position = "none") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
                                                                                                                    


firm_bac <- plot_grid(firm_plot_final, bac_plot_final, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1, rel_widths = c(2, 2))

fb_combined <- draw_plot(firm_bac, x = 0, y = 0, width = 3, height = 3)



ggsave("fig.pdf", plot=firm_bac,width = 5, height = 5, units = "in")




### Graphs for family level 

rumi_plot <- ggplot(family_all_ages, aes(x=scale_rumi, y=mgr)) + geom_point() + geom_smooth(method = "lm", se = TRUE) +  xlab(expression(paste(italic("Ruminococcaceae "), "Abundance"))) + ylab("Mass Gain Rate (g/day)")
rumi_plot_final <- rumi_plot + theme_linedraw(base_size = 22) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(breaks=c(-4,-3,-2,-1, 0, 1, 2, 3, 4, 5))


scale_muri_plot_UV <- scale(muri_plot_UV$muri_fam, center = TRUE, scale = TRUE)
muri_graph_UV <- ggplot(muri_plot_UV, aes(x = scale_muri_plot_UV, y = mgr)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red", se = TRUE) + xlab(expression(paste(italic("Muribaculaceae "), "Abundance"))) + ylab("")
muri_plot_final <- muri_graph_UV + theme_linedraw(base_size = 22) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Combine the plots

rumi_muri <- plot_grid(rumi_plot_final, muri_plot_final, 
                      labels = c("A", "B"),
                      ncol = 2, nrow = 1, rel_widths = c(2, 2))

rm_combined <- draw_plot(rumi_muri, x = 0, y = 0, width = 3, height = 3)



