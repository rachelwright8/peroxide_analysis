#install.packages("plater") # you only have to do this once per computer (it is already installed on lab Mac)

# load the packages
library(plater) # this lets you read in data in plate format
library(tidyverse) # this is for data wrangling and plotting

# read in data from plate 1-----
data <- read_plate("20220324_h2o2_data_plate1.csv", # change this to your file name!
                   well_ids_column = "position")
# look at data
# make sure everything matches your plate layout in the lab notebook
data

# remove any rows for empty wells (pvalues = NA, no sample)
data <- data %>% filter(!is.na(values))
data

# rename the columns something more descriptive
data <- data %>% rename("sample" = "values",
                        "fluorAU" = "values.2")
data

# subtract mean blank value
mean_blank <- data %>% filter(sample=="STDH") %>% 
  summarize(mean = mean(fluorAU))
mean_blank


data <- data %>% mutate(fluorBlank = fluorAU-as.numeric(mean_blank))
data

# calculate means, standard deviation, and coefficient of variance
data_summary <- data %>% select(sample,fluorBlank) %>%
  group_by(sample) %>% 
  summarise_all(funs(mean, sd), na.rm = TRUE) %>%
  mutate(cv = (sd/mean))

data_summary

# Find samples where the CV is too high
data_summary %>% filter(cv>0.2)

# these samples may need to be repeated
# standard A is out of range
data_summary <- data_summary %>% 
  filter(!sample=="STDA")


# make standard curve
# standard A was 10 µL of 3% (0.88M) in 240 µL ASW

stdA <- (0.1*0.88)/0.25
known_concs <- c(stdA^1/2, stdA^1/4, stdA^1/8, 
                 stdA^1/16, stdA^1/64, 0)
known_concs

standards <- c("STDB", "STDC", "STDD", "STDE", "STDG", "STDH")

curve_abs <- data_summary %>% 
  filter(sample %in% standards) %>%
  mutate(known_concs = known_concs) %>%
  select(known_concs, mean)
curve_abs

# plot the curve
lm_curve <- lm(mean ~ known_concs, data = curve_abs)
summary(lm_curve)
plot(mean ~ known_concs, data = curve_abs,
     ylab = "Mean Absorbance",
     xlab = "Known Concentration mM")
abline(lm_curve,col="red") 
legend("topleft", bty = "n",
       legend = paste("R2 =", format(summary(lm_curve)$r.squared,digits=3)))
# want to see R2 above 95%

# extract values for calculating unknown values
# need slope (estimate) and intercept

intercept <- as.numeric(coef(lm_curve)["(Intercept)"])
intercept
slope <- as.numeric(coef(lm_curve)["known_concs"])
slope

# make a new column with calculated concentrations by solving for the linear equation
# y = mx + b
# y is the absorbance
# m is the slope
# b is the intercept
# x is the unknown value that you are solving for

# rewrite the equation as x = (y-b)/m

# calculate concentrations of unknowns
peroxide_concentrations <- data_summary %>% 
  mutate("conc_mM" = (mean-intercept)/slope) %>%
  select(sample,conc_mM)
peroxide_concentrations

# remove the blank rows
data_summary <- data_summary %>% filter(!grepl("ASW",sample)) %>% 
  filter(!grepl("STD",sample))
data_summary <- data_summary %>% 
  mutate("treatment" = substr(sample,2,2),
         "genet" = substr(sample,3,3),
         "rep" = substr(sample,4,4),
         "plate" = rep("plate3",nrow(data_summary)))
data_summary

# read in data from plate 2-----
data2 <- read_plate("20220324_h2o2_data_plate2.csv", # change this to your file name!
                   well_ids_column = "position")
# look at data
# make sure everything matches your plate layout in the lab notebook
data2

# remove any rows for empty wells (pvalues = NA, no sample)
data2 <- data2 %>% filter(!is.na(values))
data2

# rename the columns something more descriptive
data2 <- data2 %>% rename("sample" = "values",
                        "fluorAU" = "values.2")
data2

# subtract mean blank value
mean_blank2 <- data2 %>% filter(sample=="STDH") %>% 
  summarize(mean = mean(fluorAU))
mean_blank2

data2 <- data2 %>% mutate(fluorBlank = fluorAU-as.numeric(mean_blank2))
data2

# calculate means, standard deviation, and coefficient of variance
data_summary2 <- data2 %>% select(sample,fluorBlank) %>%
  group_by(sample) %>% 
  summarise_all(funs(mean, sd), na.rm = TRUE) %>%
  mutate(cv = (sd/mean))

data_summary2

# Find samples where the CV is too high
data_summary2 %>% filter(cv>20)
# these samples may need to be repeated

# make standard curve
# standard A was 10 µL of 3% (0.88M) in 240 µL ASW

stdA <- (0.1*0.88)/0.25
known_concs <- c(stdA^1/2, stdA^1/4, stdA^1/8, 
                 stdA^1/16,stdA^1/32, stdA^1/64, 0)
known_concs

standards <- c("STDB", "STDC", "STDD", "STDE", "STDF", "STDG", "STDH")

curve_abs <- data_summary2 %>% 
  filter(sample %in% standards) %>%
  mutate(known_concs = known_concs) %>%
  select(known_concs, mean)
curve_abs

# plot the curve
lm_curve <- lm(mean ~ known_concs, data = curve_abs)
summary(lm_curve)
plot(mean ~ known_concs, data = curve_abs,
     ylab = "Mean Absorbance",
     xlab = "Known Concentration mM")
abline(lm_curve,col="red") 
legend("topleft", bty = "n",
       legend = paste("R2 =", format(summary(lm_curve)$r.squared,digits=3)))
# want to see R2 above 95%

# extract values for calculating unknown values
# need slope (estimate) and intercept

intercept <- as.numeric(coef(lm_curve)["(Intercept)"])
intercept
slope <- as.numeric(coef(lm_curve)["known_concs"])
slope

# make a new column with calculated concentrations by solving for the linear equation
# y = mx + b
# y is the absorbance
# m is the slope
# b is the intercept
# x is the unknown value that you are solving for

# rewrite the equation as x = (y-b)/m

# calculate concentrations of unknowns
peroxide_concentrations2 <- data_summary2 %>% 
  mutate("conc_mM" = (mean-intercept)/slope) %>%
  select(sample,conc_mM)
peroxide_concentrations2

# remove the blank rows
data_summary2 <- data_summary2 %>% 
  filter(!grepl("ASW",sample)) %>% 
  filter(!grepl("STD",sample))

data_summary2 <- data_summary2 %>% 
  mutate("treatment" = substr(sample,2,2),
         "genet" = substr(sample,3,3),
         "rep" = substr(sample,4,4),
         "plate" = rep("plate2",nrow(data_summary2)))
data_summary2



# read in data from plate 3-----
data3 <- read_plate("20220324_h2o2_data_plate3.csv", # change this to your file name!
                    well_ids_column = "position")
# look at data
# make sure everything matches your plate layout in the lab notebook
data3

# remove any rows for empty wells (pvalues = NA, no sample)
data3 <- data3 %>% filter(!is.na(values))
data3

# rename the columns something more descriptive
data3 <- data3 %>% rename("sample" = "values",
                          "fluorAU" = "values.2")
data3

# subtract mean blank value
mean_blank3 <- data3 %>% filter(sample=="STDH") %>% 
  summarize(mean = mean(fluorAU))
mean_blank3

data3 <- data3 %>% 
  mutate(fluorBlank = fluorAU-as.numeric(mean_blank3))
data3

# calculate means, standard deviation, and coefficient of variance
data_summary3 <- data3 %>% select(sample,fluorBlank) %>%
  group_by(sample) %>% 
  summarise_all(funs(mean, sd), na.rm = TRUE) %>%
  mutate(cv = (sd/mean))

# Find samples where the CV is too high
data_summary3 %>% filter(cv>20)
# these samples may need to be repeated

# make standard curve
# standard A was 10 µL of 3% (0.88M) in 240 µL ASW

stdA <- (0.1*0.88)/0.25
known_concs <- c(stdA^1/2, stdA^1/4, stdA^1/8, 
                 stdA^1/16,stdA^1/32, stdA^1/64, 0)
known_concs

standards <- c("STDB", "STDC", "STDD", "STDE", "STDF", "STDG", "STDH")

curve_abs <- data_summary3 %>% 
  filter(sample %in% standards) %>%
  mutate(known_concs = known_concs) %>%
  select(known_concs, mean)
curve_abs

# plot the curve
lm_curve <- lm(mean ~ known_concs, data = curve_abs)
summary(lm_curve)
plot(mean ~ known_concs, data = curve_abs,
     ylab = "Mean Absorbance",
     xlab = "Known Concentration mM")
abline(lm_curve,col="red") 
legend("topleft", bty = "n",
       legend = paste("R2 =", format(summary(lm_curve)$r.squared,digits=3)))
# want to see R2 above 95%

# extract values for calculating unknown values
# need slope (estimate) and intercept

intercept <- as.numeric(coef(lm_curve)["(Intercept)"])
intercept
slope <- as.numeric(coef(lm_curve)["known_concs"])
slope

# make a new column with calculated concentrations by solving for the linear equation
# y = mx + b
# y is the absorbance
# m is the slope
# b is the intercept
# x is the unknown value that you are solving for

# rewrite the equation as x = (y-b)/m

# calculate concentrations of unknowns
peroxide_concentrations3 <- data_summary3 %>% 
  mutate("conc_mM" = (mean-intercept)/slope) %>%
  select(sample,conc_mM)
peroxide_concentrations3

# remove the blank rows
data_summary3 <- data_summary3 %>% 
  filter(!grepl("ASW",sample)) %>% 
  filter(!grepl("STD",sample))

data_summary3 <- data_summary3 %>% 
  mutate("treatment" = substr(sample,2,2),
         "genet" = substr(sample,3,3),
         "rep" = substr(sample,4,4),
         "plate" = rep("plate3",nrow(data_summary3)))
data_summary3

# combined all plates

all_data <- rbind(data_summary,
                      data_summary2,
                      data_summary3)

all_data_anemones <- all_data %>% filter(!genet=="A")
all_data_algae <- all_data %>% filter(genet=="A")

# across genets

all_data %>% filter(genet=="B" & treatment=="N")

all_data  %>%
  ggplot(aes(x=genet,y=mean,fill=treatment))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_manual(values=c("lightblue","darkorange","grey"))+
  scale_x_discrete(labels = levels(as.factor(all_data$genet))) +
  ylab("Peroxide Concentration (mM)")+
  xlab("Genet")+
  ggtitle("Peroxide Concentrations")

# across treatments
all_data_anemones  %>%
  ggplot(aes(x=treatment,y=mean,fill=treatment))+
  geom_violin()+
  theme_classic()+
  theme(legend.position="none")+
  # scale_fill_manual(values=c("lightblue","darkorange","grey"))+
  # scale_x_discrete(labels = c("Cold","Hot","Control")) +
  scale_fill_manual(values=c("lightblue","grey"))+
  scale_x_discrete(labels = c("Cold","Control")) +
  ylab("Peroxide Concentration (mM)")+
  xlab("Stress Treatment")+
  ggtitle("Peroxide in Anemones")

all_data_anemones  %>%
  ggplot(aes(x=treatment,y=mean,fill=treatment))+
  geom_boxplot()+
  geom_point()+
  theme_classic()+
  theme(legend.position="none")+
  # scale_fill_manual(values=c("lightblue","darkorange","grey"))+
  # scale_x_discrete(labels = c("Cold","Hot","Control")) +
  scale_fill_manual(values=c("lightblue","grey"))+
  scale_x_discrete(labels = c("Cold","Control")) +
  ylab("Peroxide Concentration (mM)")+
  xlab("Stress Treatment")+
  ggtitle("Peroxide in Anemones")

all_data_algae  %>%
  ggplot(aes(x=treatment,y=mean,fill=treatment))+
  geom_boxplot()+
  geom_point()+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("lightblue","darkorange","grey"))+
  scale_x_discrete(
    labels = c("Cold","Hot","Control")) +
  ylab("Peroxide Concentration (mM)")+
  xlab("Stress Treatment")+
  ggtitle("Peroxide in Algal Culture")

# stats ---

# algae cultures stats

all_data_algae
all_data_algae$treatment <- fct_relevel(all_data_algae$treatment,"N")
lmCulture <- lm(mean ~ treatment, data=all_data_algae)
summary(lmCulture)
anova(lmCulture)

aovCulture <- aov(mean ~ treatment, data=all_data_algae)
summary(aovCulture)
#               Df   Sum Sq   Mean Sq F value  Pr(>F)   
# treatment    2 0.005259 0.0026297   10.49 0.00232 **
# Residuals   12 0.003007 0.0002506 

TukeyHSD(x=aovCulture, conf.level=0.95)
#     diff        lwr       upr     p adj
# C-N -0.01666667 -0.041050717 0.007717383 0.2036701
# H-N  0.03461111  0.004746871 0.064475351 0.0235131
# H-C  0.05127778  0.021413538 0.081142018 0.0016886

# anemones cultures stats

all_data_anemones_nohot
all_data_anemones_nohot$treatment <- fct_relevel(all_data_anemones_nohot$treatment,"N")
lmAnemone <- lm(mean ~ treatment, data=all_data_anemones_nohot)
summary(lmAnemone)
anova(lmAnemone)

aovAnemone <- aov(mean ~ treatment, data=all_data_anemones_nohot)
summary(aovAnemone)
#             Df  Sum Sq   Mean Sq F value Pr(>F)  
# treatment    1 0.00208 0.0020803   6.576 0.0137 *
# Residuals   46 0.01455 0.0003164  

TukeyHSD(x=aovAnemone, conf.level=0.95)
#         diff        lwr        upr     p adj
# C-N 0.01316667 0.00283149 0.02350184 0.0136667

# save data -----
all_data_20220324 <- all_data
save(all_data_20220324, file = "alldata_20220324.Rdata")
