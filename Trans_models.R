#----Load libraries and check data----
library(ggpubr)
library(PerformanceAnalytics)
library(tidyverse)
library(stats)
library(MuMIn)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(ggplot2)
library(broom.mixed)
library(car)
library(Performance)

# Read in script
Trans_2 <-read.csv('Transloc_sum_V5.csv', header = TRUE)
summary(Trans_2)

# Standardize key variables
Trans_2 <- Trans_2 %>%
  mutate(Transmitters_std = as.numeric(scale(Tranmitters))) %>%
  mutate(Training_std = as.numeric(scale(Training))) %>%
  mutate(Conspecifics_std = as.numeric(scale(Conspecifics))) %>%
  mutate(Total_release_std = as.numeric(scale(Total_release)))

# Check correlations
correlation_subset <- Trans_2[, c(10,12,14,15,30,31,32,33)]
chart.Correlation(correlation_subset, histogram=TRUE)

#----Run models on translocation survival----
f1 <- glm(cbind(Transloc, Total_r - Transloc) ~ Total_release_std, family="binomial", data=Trans_2)
f2 <- glm(cbind(Transloc, Total_r - Transloc) ~ Transmitters_std, family="binomial", data=Trans_2)
f3 <- glm(cbind(Transloc, Total_r - Transloc) ~ Training_std, family="binomial", data=Trans_2)
f4 <- glm(cbind(Transloc, Total_r - Transloc) ~ Total_release_std + Transmitters_std + Training_std , family="binomial", data=Trans_2)
f5 <- glm(cbind(Transloc, Total_r - Transloc) ~ 1 , family="binomial", data=Trans_2)

summary(f4)
summary(f1)
vif(f4)

# Run model selection
fmList1<-model.sel(f1=f1, f2=f2, f3=f3, f4=f4, f5=f5)
fmList1

# Check residuals
sim_f4 <- simulateResiduals(fittedModel = f4, n = 1000)
plot(sim_f4)
testDispersion(sim_f4)
testUniformity(sim_f4)
testOutliers(sim_f4)

which(sim_f4$outliers)
plot(sim_f4)

# Make model predictions plot
newdata <- data.frame(
  Total_release_std = seq(min(Trans_2$Total_release_std),
                          max(Trans_2$Total_release_std),
                          length = 100),
  Transmitters_std = mean(Trans_2$Transmitters_std, na.rm = TRUE),
  Training_std = mean(Trans_2$Training_std, na.rm = TRUE))

pred <- predict(f4, newdata = newdata, type = "response", se.fit = TRUE)

newdata$fit <- pred$fit
newdata$se <- pred$se.fit
newdata$lwr <- newdata$fit - 1.96 * newdata$se
newdata$upr <- newdata$fit + 1.96 * newdata$se
Trans_2$prop_survival <- Trans_2$Transloc / Trans_2$Total_r

f4plot <- ggplot() +
  geom_point(data = Trans_2,
             aes(x = Total_release_std, y = prop_survival),
             size = 2, alpha = 0.7) +
  geom_line(data = newdata,
            aes(x = Total_release_std, y = fit),
            size = 1.2, colour = "black") +
  geom_ribbon(data = newdata,
              aes(x = Total_release_std, ymin = lwr, ymax = upr),
              alpha = 0.2, fill = "black") +
  labs(x = "Release group size (standardized)",
       y = "Translocation survival probability") +
  theme_classic()
f4plot

#----Run models on season survival----
b1 <- glm(cbind(Season_survive, Total_r - Season_survive) ~ Total_release_std, family="binomial", data=Trans_2)
b2 <- glm(cbind(Season_survive, Total_r - Season_survive) ~ Transmitters_std, family="binomial", data=Trans_2)
b3 <- glm(cbind(Season_survive, Total_r - Season_survive) ~ Training_std, family="binomial", data=Trans_2)
b4 <- glm(cbind(Season_survive, Total_r - Season_survive) ~ Total_release_std + Transmitters_std + Training_std , family="binomial", data=Trans_2)
b5 <- glm(cbind(Season_survive, Total_r - Season_survive) ~ 1 , family="binomial", data=Trans_2)

summary(b2)
summary(b4)
vif(b4)

# Run model selection
fmList<-model.sel(b1=b1, b2=b2, b3=b3, b4=b4, b5=b5)
fmList
r2(b1)

# Check residuals
sim_b2 <- simulateResiduals(fittedModel = b2, n = 1000)
plot(sim_b2)
testDispersion(sim_b2)
testUniformity(sim_b2)
testOutliers(sim_b2)

which(sim_b2$outliers)
plot(sim_b2)

# Make model predictions plot
newdata <- data.frame(
  Transmitters_std = seq(min(Trans_2$Transmitters_std),
                         max(Trans_2$Transmitters_std),
                         length = 100))

pred <- predict(b2, newdata = newdata, type = "response", se.fit = TRUE)

newdata$fit <- pred$fit
newdata$se <- pred$se.fit
newdata$lwr <- newdata$fit - 1.96 * newdata$se
newdata$upr <- newdata$fit + 1.96 * newdata$se
Trans_2$prop_season_survival <- Trans_2$Season_survive / Trans_2$Total_r
b2plot <- ggplot() +
  geom_point(data = Trans_2,
             aes(x = Transmitters_std, y = prop_season_survival),
             size = 2, alpha = 0.7) +
  geom_line(data = newdata,
            aes(x = Transmitters_std, y = fit),
            colour = "black", linewidth = 1.2) +
  geom_ribbon(data = newdata,
              aes(x = Transmitters_std, ymin = lwr, ymax = upr),
              fill = "black", alpha = 0.2) +
  labs(x = "Transmitters (standardized)",
       y = "Predicted seasonal survival probability") +
  theme_classic()
b2plot

# put model predictions plots together
library(patchwork)

f4plot + b2plot +
  plot_annotation(tag_levels = "a")

#----Plot coefficent estimates for both model sets----
## Extract coefficient estimates
f4_coefs <- tidy(f4, effects = "fixed") %>%
  mutate(model = "f4")
b2_coefs <- tidy(b2, effects = "fixed") %>%
  mutate(model = "b2")
coef_data <- bind_rows(f4_coefs, b2_coefs)

coef_data <- coef_data %>%
  rename(Variable = term) %>%
  rename(Coefficent_estimate = estimate) %>%
  mutate(model = dplyr::recode(model, "f4" = "Translocation_survival")) %>%
  mutate(model = dplyr::recode(model, "b2" = "Season_survival"))


p <- ggplot(coef_data, aes(x = Coefficent_estimate, y = Variable)) +
  geom_point() +
  geom_errorbarh(aes(xmin = Coefficent_estimate - 1.96*std.error,
                     xmax = Coefficent_estimate + 1.96*std.error),
                 height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~model) +
  theme_bw()
p

#----Make survival plots----
# Make plot of proportion that survive the season and poroportion that survive the translocation
library(ggplot2)
p = ggplot(Trans_2,aes(x=Prop_Season_survival, y=Prop_Transloc))
p = p + geom_point()
p = p + theme(panel.background = element_rect(fill = 'white', colour = 'darkgrey'))
p = p + geom_smooth(method=lm, level=0.99, col='black')
p = p + theme(panel.grid.major = element_line(colour = 'white'))
p = p + labs(x= "Proportion of birds that survive the translocation", y= "Proportion of birds that survive the season")
p = p + theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))
p = p + theme(text = element_text(size = 20))    
print(p)

# Make plot of proportion that survive the season and poroportion that return the following season

p = ggplot(Trans_2,aes(x=Prop_Returns,y=Prop_Season_survival))
p = p + geom_smooth(method=lm, level=0.99, col='black')
p = p + geom_point()
p = p + theme(panel.background = element_rect(fill = 'white', colour = 'darkgrey'))
p = p + theme(panel.grid.major = element_line(colour = 'white'))
p = p + labs(x= "Proportion of birds that return the following season", y= "Proportion of birds that survive the season")
p = p + theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))
p = p + theme(text = element_text(size = 20))  
print(p)

# Translocation survival and total release size
p = ggplot(Trans_2,aes(x=Transloc,y=Total_release))
p = p + geom_point()
p = p + geom_smooth(method=lm, level=0.99, col='black')
p = p + theme(panel.background = element_rect(fill = 'white', colour = 'darkgrey'))
p = p + theme(panel.grid.major = element_line(colour = 'white'))
p = p + labs(x= "Number of birds surviving the translocation", y= "Total number of released birds")
p = p + theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))
print(p)

# Release size and proportion that return
p = ggplot(comp,aes(x=Prop_Returns,y=Total_release))
p = p + geom_smooth(method=lm, level=0.99, col='black')
p = p + geom_point()
p = p + theme(panel.background = element_rect(fill = 'white', colour = 'darkgrey'))
p = p + theme(panel.grid.major = element_line(colour = 'white'))
p = p + labs(x= "Proportion of birds that return the following season", y= "Total release size")
p = p + theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))
print(p)

## Plot survival differences between locations
Ave <-read.csv('Average_success2.csv', header = TRUE)
head(Ave)
p = ggplot(data=Ave, aes(x=Location, y=Percentage, fill = Value)) + geom_bar(stat="identity", position="dodge")
p = p + scale_color_manual(values=c("#333333", "#999999", "#CCCCCC")) + scale_fill_manual(values=c("#333333","#999999", "#CCCCCC"))
p = p + geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.2,position=position_dodge(.9))
p = p + theme(panel.background = element_rect(fill = 'white', colour = 'white'))
p = p + theme(panel.border = element_blank(), axis.line = element_line())
p = p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p = p + labs(x= "Location", y= "Average percentage")
p = p + theme(axis.text.x = element_text(colour = "black", size=15, margin= margin(t=0,r=20,b=0,l=0)), 
              axis.text.y = element_text(colour = "black", size=15, margin = margin(t=20,r=0,b=20,l=0)))
p = p + theme(axis.title=element_text(size=15))
p = p + theme(legend.text=element_text(size=14))
p = p + theme(legend.title=element_text(size=14))
print(p)

#----Plot consistency between data sources----
Ave <-read.csv('Consistency.csv', header = TRUE)
head(Ave)
p = ggplot(data=Ave, aes(x=Year, y=Birds_released, fill = Source)) + geom_bar(stat="identity", position="dodge")
#p = p + scale_color_manual(values=c("#333333", "#999999", "#CCCCCC")) + scale_fill_manual(values=c("#333333","#999999", "#CCCCCC"))
#p = p + geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.2,position=position_dodge(.9))
p = p + theme(panel.background = element_rect(fill = 'white', colour = 'white'))
p = p + theme(panel.border = element_blank(), axis.line = element_line())
#p = p + scale_x_continuous(breaks = 0:11,
#labels = paste0(c("1994", "1996", "1999", "2000", "2001", "2002", "2003"), "2004", "2005", "2007", "2008", "2009"), "year")
#p = p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p = p + labs(x= "Year", y= "Number of birds released")
p = p + theme(axis.text.x = element_text(colour = "black", size=15, margin= margin(t=0,r=20,b=0,l=0)), 
              axis.text.y = element_text(colour = "black", size=15, margin = margin(t=20,r=0,b=20,l=0)))
p = p + theme(axis.title=element_text(size=15))
p = p + scale_fill_brewer(palette="YlGnBu")
p = p + theme(legend.text=element_text(size=14))
p = p + theme(legend.title=element_text(size=14))
print(p)
