rm(list = ls())
#packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, ggplot2, tidyr, readxl, tidyverse, asreml, ggrepel) 

#Get data
df <- read.csv(file = "~/Documents/INIA/analysis_ainhoa/datafiles/CUR_data.csv", row.names=NULL)
str(df)
names(df)

df$Height <- as.numeric(df$Height)

str(df)
head(df); tail(df)
#*********************************
# Review data
mean_long <- df %>%
  group_by(PROV, REP, YP) %>%
  summarise(mean_height = mean(Height, na.rm = TRUE), .groups = "drop")

# as.factor() required for ASRemlR
df$COL <- as.factor(df$COL)
df$ROW <- as.factor(df$ROW)
df$PROV <- as.factor(df$PROV)
df$REP <- as.factor(df$REP)
df$REG <- as.factor(df$REG)
df$POS <- as.factor(df$POS)
str(df)
#df$Age <- as.factor(df$Age)

# Data evaluation prior to analysis
table(df$PROV, df$REP)
table(df$PROV, df$POS)
hist(df$Height)
boxplot(Height ~ Age, data = df)
boxplot(Height ~ YP, data = df)
meanH_A <- aggregate(Height ~ Age, mean, data = df)
meanH_AP <- aggregate(Height ~ PROV + Age, FUN = function(x) mean(x, na.rm = TRUE), data = df)

labels_data <- meanH_AP %>%
  group_by(PROV) %>%
  filter(Age == max(Age))  # Select the last available Age for each PROV

H_AP_plot <- ggplot(data=meanH_AP, aes(x=Age, y=Height, group = PROV)) +
  geom_line(linewidth=0.7, aes(color=PROV)) + 
  #geom_text(data=labels_data, aes(label=PROV), hjust=-0.1, size=4) + 
  geom_text_repel(data=labels_data, aes(label=PROV), size=4, direction = "y", nudge_x = 2) +  
  ylab("Height (cm)") + xlab("Age (years)") 

H_AP_plot

#********************************************
# Check impact of different beta*(1/(x+c)) in linearity
df_imp <- data.frame(Age = 5:11)
df_imp$iAge_0 <- 1 / (df_imp$Age)
df_imp$iAge_3 <- 1 / (df_imp$Age + 3)
df_imp$iAge_5 <- 1 / (df_imp$Age + 5)
df_imp$iAge_7 <- 1 / (df_imp$Age + 7)
df_imp$iAge_10 <- 1 / (df_imp$Age + 10)
df_imp$iAge_20 <- 1 / (df_imp$Age + 20)

df_imp_long <- df_imp %>%
  pivot_longer(cols = starts_with("iAge_"), names_to = "Transformation", values_to = "iAge_Value")

ggplot(df_imp_long, aes(x = Age, y = iAge_Value, color = Transformation)) +
  geom_line(size = 1) +
  labs(
    title = "Impact of Different Transformations on iAge",
    x = "Age (years)",
    y = "Transformed Age (iAge)"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "bottom")

#************* Data Linearisation ************************

df$lnY <- log(df$Height+0)
df$iAge <- 1/(df$Age+5)  
hist(df$lnY)
boxplot(lnY~iAge, data=df)

mean_lnY_AP <- aggregate(lnY ~ PROV + iAge, FUN = function(x) mean(x, na.rm = TRUE), data = df)

lnY_AP_plot <- ggplot(data=mean_lnY_AP, aes(x=iAge, y=lnY, group = PROV)) +
  geom_line(linewidth=0.7, aes(color=PROV)) + 
  ylab("lnY") + xlab("iAge") 

lnY_AP_plot

#**************** RR models ******************
# as.factor() required for ASRemlR
df$COL <- as.factor(df$COL)
df$ROW <- as.factor(df$ROW)
df$PROV <- as.factor(df$PROV)
df$REP <- as.factor(df$REP)
df$REG <- as.factor(df$REG)
df$POS <- as.factor(df$POS)
df$YP <- as.factor(df$YP)
#df$Age <- as.factor(df$Age)

table(df$PROV, df$REP)

#****** Model 1, independent RE *****
rr_mod1 <- asreml(fixed = lnY ~ 1 + REP + iAge,
              random = ~ PROV + PROV:iAge,
              residual = ~ idv(units),
              data = df)
wald(rr_mod1)
summary(rr_mod1)$varcomp
rr_mod1$loglik 
summary(rr_mod1)$aic
summary(rr_mod1)$bic

#****** Model 2, correlated RE, independent res *****
length(levels(df$PROV))  ## 19 levels
length(levels(df$REG))   
asreml.options(maxit=25, gammaPar = FALSE, ai.sing = TRUE, fail="soft", workspace = "4gb")
rr_mod2 <- asreml(fixed = lnY ~ 1 + REP + iAge,
              random = ~str(~PROV + PROV:iAge, ~corgh(2):id(16)),
              residual = ~idv(units),
              data = df)

wald(rr_mod2)
summary(rr_mod2)$varcomp                   
rr_mod2$loglik  
summary(rr_mod2)$aic
summary(rr_mod2)$bic

#model comparisson
lrt.asreml(rr_mod1,rr_mod2,boundary = FALSE)

#****** Model 3, correlated RE and correlated res *****
length(levels(df$REG))
length(levels(df$YP))
df <- df[order(df$REG, decreasing = FALSE),]
df <- df[order(df$YP, decreasing = FALSE),]
asreml.options(maxit=25, gammaPar = FALSE, ai.sing = TRUE, fail="soft", workspace = "4gb")

initR <- c(0.8, 0.09)
rr_mod3 <- asreml(fixed = lnY ~ 1 + REP + iAge,
              random = ~str(~PROV + PROV:iAge, ~corgh(2):id(16)),
              residual = ~ar1v(YP, init=initR):REG,
              data = df)

wald(rr_mod3)
summary(rr_mod3)$varcomp
rr_mod3$loglik 
summary(rr_mod3)$aic
summary(rr_mod3)$bic
lrt.asreml(rr_mod2,rr_mod3,boundary = FALSE)
plot(rr_mod3)


#BLUE
BLUE.mod3 <- summary(rr_mod3, coef = TRUE)$coef.fixed

# Getting alpha -> overall intercept
#1/(0+5) 
alpha_int <- predict.asreml(rr_mod3, classify='(Intercept):iAge', levels=list(iAge=0))$pvals  #overall intercept
alpha_int  # 8.466253

# Getting beta -> overall slope
BLUE.mod3
beta_slope <- BLUE.mod3[6,1]  

#BLUPS
BLUP.mod3 <- summary(rr3, coef = TRUE)$coef.random
head(BLUP.mod3); tail(BLUP.mod3)
dim(BLUP.mod3)

BLUP_mod3a <- as.data.frame(BLUP.mod3[1:18,])  #PROV
BLUP_mod3b <- as.data.frame(BLUP.mod3[19:36,])  #PROV:iAge

RRcoef <- data.frame(PROV = rownames(BLUP_mod3a), BLUP_mod3a = BLUP_mod3a$solution, BLUP_mod3b = BLUP_mod3b$solution)

# Getting a and b per PROV
RRcoef$a_int <- alpha_int$predicted.value + RRcoef$BLUP_mod3a
RRcoef$b_slope <- beta_slope + RRcoef$BLUP_mod3b
head(RRcoef); tail(RRcoef)

# Predictions for age 40 (rotation age)
iAge_40 <- 1/(40+5)  #iAge = 1/(Age +5)
pred_40Age <- predict.asreml(rr3, classify = "PROV:iAge", levels = list(iAge=iAge_40))$pvals

df_coef <- RRcoef  #To not modify RRcoef 
df_coef$pred_Age40 <- exp(pred_40Age$predicted.value)-0  # extrapolated and back-transformed


# Predictions for ages
m_iAges <- unique(df$iAge)
m_iAges <- 1/(seq(7,26,1)+5)
pred_iAges <- predict.asreml(rr_mod3, classify = "PROV:iAge", levels = list(iAge=m_iAges))$pvals
pred_Ages <- as.data.frame(pred_iAges)

# Back transforming 
pred_Ages$Age <- round(1/pred_Ages$iAge - 5.0)
pred_Ages$pred_Height <- exp(pred_Ages$predicted.value)
head(pred_Ages); tail(pred_Ages)


