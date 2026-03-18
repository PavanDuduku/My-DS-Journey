#----------Section 1-------------------
# Set working directory
setwd(dirname(file.choose()))
getwd()

# Get data from CSV file
covid.deaths <- read.csv("cw-data.csv", stringsAsFactors = FALSE)
#row.names(covid.deaths) <- covid.deaths[, "LA_name"]
#covid.deaths <- covid.deaths[,-1]
head(covid.deaths)
str(covid.deaths)

# Check for missing data
library(Amelia)
missmap(covid.deaths, col = c("black", "grey"), legend = FALSE, x.cex = 0.6, margins = c(6, 4))


#----------Section 2-------------------
# Check dependent variable's distribution

boxplot(covid.deaths$Total_Deaths_1000, xlab="Total Covid deaths per 1000", ylab="Count")

# Label outliers for each outlier in boxdata
boxdata <- boxplot(covid.deaths$Total_Deaths_1000, xlab="Total Covid deaths per 1000", ylab="Count")
text_pos_flag = TRUE
for (i in 1:length(boxdata$group)) {
  text (boxdata$group[i], boxdata$out[i], 
        labels = sprintf("%s (%.2f)", covid.deaths$LA_name[covid.deaths$Total_Deaths_1000 == boxdata$out[i]], boxdata$out[i]), 
        pos = ifelse(text_pos_flag, 4, 2), cex = 0.7)
  text_pos_flag = !text_pos_flag
}
rm(boxdata)

# Isles of Scilly is an outlier as it has 0 deaths, we'll remove it & other columns for further analysis
covid.deaths <- covid.deaths[covid.deaths$LA_name != "Isles of Scilly", !(names(covid.deaths) %in% c("LA_name", "LA_code", "Total_Population"))]
attach(covid.deaths)

# Frequency histogram
hist(Total_Deaths_1000, freq = TRUE)

# Probability density histogram 
hist(Total_Deaths_1000, freq = FALSE, ylim = c(0,0.8))
# Add a density curve
lines(density(sort(Total_Deaths_1000)))

# Add a Normal curve
xfit <- seq(from = min(Total_Deaths_1000), to = max(Total_Deaths_1000), by = 0.01)
yfit = dnorm(xfit, mean(Total_Deaths_1000), sd(Total_Deaths_1000))
lines(xfit, yfit, lty = "dotted")
rm(xfit, yfit)

# Add a rug plot
rug(Total_Deaths_1000)

# Add a legend
legend("topright", legend = c("Density curve", "Normal curve"),
       lty = c("solid", "dotted"), cex = 0.7)

# Test dependent variable for normality graphically
qqnorm(Total_Deaths_1000)
qqline(Total_Deaths_1000, col = 2) ## red color

# Significance testing for normality
ks.test(Total_Deaths_1000,"pnorm", mean(Total_Deaths_1000), sd(Total_Deaths_1000)) #0.9339
shapiro.test(Total_Deaths_1000) #0.2675


#----------Section 3-------------------
# Check independent variables for normality (to decide which significance test)

independent_vars <- setdiff(names(covid.deaths), "Total_Deaths_1000")
normality.df <- data.frame( "Variable" = character(),
                            "K-S_p_value" = numeric(),
                            "ShapiroWilk_p_value" = numeric())
for (var in independent_vars) {
  x <- as.numeric(covid.deaths[,var])
  ks.res <- ks.test(x ,"pnorm", mean(x), sd(x))
  ks.res.p <- ifelse(ks.res$p.value < 1e-4,
                     format(ks.res$p.value, scientific = TRUE),
                     format(ks.res$p.value, scientific = FALSE))
  sh_w.res <- shapiro.test(x)
  sh_w.res.p <- ifelse(sh_w.res$p.value < 1e-4,
                       format(sh_w.res$p.value, scientific = TRUE),
                       format(sh_w.res$p.value, scientific = FALSE))
  normality.df[nrow(normality.df) +1, ] <- list(var, ks.res.p, sh_w.res.p)
}

# In normality-check functions like Shapiro Wilk & K-S test, null hypothesis = The data do not deviate from normality => follows normal distribution
# We'll do a visual check on normality ambiguous variables (either of the two p-values > 0.05)
ambiguos_vars = c("X15_and_under", "X35_to_49", "X50_to_74", 
                  "Very_good_health", "Good_health", "Bad_health", "Very_bad_health",
                  "HH_no_depr", "HH_Depr_2_dim")
library(psych)
for (var in ambiguos_vars) {
  x <- as.numeric(covid.deaths[,var])
  qqnorm(x, main = var)
  qqline(x, col = 2)
  readline("Press Enter for next plot...") # Needs interaction with console
  # Probability density histogram
  hist(x, freq = FALSE, main = paste("Histogram of" , var))
  # Add a density curve
  lines(density(sort(x)))
  # Add a Normal curve
  xfit <- seq(from = min(x), to = max(x), by = 0.01)
  yfit = dnorm(xfit, mean(x), sd(x))
  lines(xfit, yfit, lty = "dotted")
  rm(xfit, yfit)
  summary_stats <- describe(x)
  print(summary_stats[, c("skew", "kurtosis")])
  readline("Press Enter for next variable's plot...") # Needs interaction with console
}

# Upon mathematically & visually examining, below are the reasonably normally distributed variables
normal_dist.vars = c("X35_to_49", "X50_to_74" ,"X75_and_over",
                     "Very_good_health", "Fair_health", "Bad_health", "Very_bad_health",
                     "HH_no_depr", "HH_Depr_1_dim", "HH_Depr_2_dim")
non_normal_dist.vars = setdiff(names(covid.deaths), c(normal_dist.vars, "Total_Deaths_1000"))


#----------Section 4-------------------
# Significance test between independent variables & dependent variable 

# Null hypothesis = There is no significant difference between chosen independent variable & Total_Deaths_1000
#                   (OR) There's no association between dependent & Independent variable (OR) They are statistically independent
significance.df <- data.frame("Variable" = character(),
                              "first_test_p_value" = numeric(),
                              "second_test_p_value" = numeric(),
                              "effect_size" = numeric())
library(rcompanion)
for (var in non_normal_dist.vars) {
  mytable <- as.data.frame(cbind(covid.deaths[,var], covid.deaths[,"Total_Deaths_1000"]))
  chisq.res <- chisq.test(mytable)
  chisq.res.p <- ifelse(chisq.res$p.value < 1e-4,
                        format(chisq.res$p.value, scientific = TRUE),
                        format(chisq.res$p.value, scientific = FALSE))
  cramer.res <- cramerV(as.matrix(mytable), conf = 0.95)
  significance.df[nrow(significance.df) +1, ] <- list(var, chisq.res.p, NA, cramer.res)
}

library(lsr)
for (var in normal_dist.vars) {
  x <- as.numeric(covid.deaths[,var])
  t.test.res <- t.test(x, Total_Deaths_1000, paired=TRUE, var.equal = FALSE)
  t.test.res.p <- ifelse(t.test.res$p.value < 1e-4,
                         format(t.test.res$p.value, scientific = TRUE),
                         format(t.test.res$p.value, scientific = FALSE))
  wil.test.res <- wilcox.test(x, Total_Deaths_1000, paired = TRUE)
  wil.test.res.p <- ifelse(wil.test.res$p.value < 1e-4,
                           format(wil.test.res$p.value, scientific = TRUE),
                           format(wil.test.res$p.value, scientific = FALSE))
  cohens.res <- cohensD(x, Total_Deaths_1000, method = "unequal")
  significance.df[nrow(significance.df) +1, ] <- list(var, t.test.res.p, wil.test.res.p, cohens.res)
}

# Statistically significant (=> p<0.05) variables are 
#   (X35_to_49, X50_to_74, X75_and_over,
#    Very_good_health, Fair_health, Bad_health, Very_bad_health
#    Asian_All, Black_All, Other,
#    HH_no_depr, HH_Depr_1_dim, HH_Depr_2_dim, 
#    x_gt_1p5_to_2)


#----------Section 5---------
# Check correlation with dependent variable 

#library(psych)
#cor.ci(covid.deaths, method = "spearman")

#library(corrgram)
##corrgram works best with Pearson correlation
#corrgram(covid.deaths, order=FALSE, cor.method = "pearson", lower.panel=panel.conf,
#         upper.panel=panel.pie, text.panel=panel.txt, main="variables")

library(corrplot)
# Correlation matrix for normal variables
cor.matrix1 <- cor(covid.deaths[,c(normal_dist.vars, "Total_Deaths_1000")], use = "pairwise.complete.obs", method = "pearson")
corrplot(cor.matrix1, type = "upper", tl.col = "black", tl.srt = 45)
# Correlation matrix for non-normal variables
cor.matrix2 <- cor(covid.deaths[,c(non_normal_dist.vars, "Total_Deaths_1000")], use = "pairwise.complete.obs", method = "spearman") #rank based
corrplot(cor.matrix2, type = "upper", tl.col = "black", tl.srt = 45)

## We'll take a detour & cleanup some variables in people per room
#detach(covid.deaths)
#covid.deaths <- within(covid.deaths, x_gt_eq_1 <- x_eq_1 + x_gt_1_to_1p5 + x_gt_1p5_to_2 + x_gt_2)
#covid.deaths <- covid.deaths[,setdiff(names(covid.deaths), c("x_eq_1", "x_gt_1_to_1p5", "x_gt_1p5_to_2", "x_gt_2", "Total_Deaths_1000"))]
#covid.deaths <- data.frame(covid.deaths, Total_Deaths_1000)
#attach(covid.deaths)
#cor.matrix <- cor(covid.deaths, use = "pairwise.complete.obs", method = "spearman")
#corrplot(cor.matrix, type = "upper", tl.col = "black", tl.srt = 45)
#cor.matrix2 <- cor(covid.deaths, use = "pairwise.complete.obs", method = "pearson")
#corrplot(cor.matrix2, type = "upper", tl.col = "black", tl.srt = 45)


# General health & Deprivation themes stood out in correlation with Covid death count
# while Ethnicity & people per room in household doesn't show any correlation.
# Age shows very faint correlation

cor.matrix1[,ncol(cor.matrix1)]
cor.matrix2[,ncol(cor.matrix2)]

# We'll take correlation coefficient threshold as 0.2, so resulting variables from correlation analysis are
#   (X15_and_under,
#    Very_good_health, Fair_health, Bad_health, Very_bad_health,
#    HH_no_depr, HH_Depr_1_dim, HH_Depr_2_dim, HH_Depr_3_dim)


#----------Section 6--------
# We'll formulate a dataset that is union of results from Significance tests & Correlation 
covid.deaths.df <- covid.deaths[,c("X15_and_under", "X35_to_49", "X50_to_74", "X75_and_over",
                                   "Very_good_health", "Fair_health", "Bad_health", "Very_bad_health",
                                   "Asian_All", "Black_All", "Other",
                                   "HH_no_depr", "HH_Depr_1_dim", "HH_Depr_2_dim", "HH_Depr_3_dim",
                                   "x_gt_1p5_to_2",
                                   "Total_Deaths_1000")]

# Exclude dependent variable 
covid.deaths.df1 <- covid.deaths.df[, setdiff(names(covid.deaths.df), "Total_Deaths_1000")]

model1 <- lm(Total_Deaths_1000 ~ ., data = covid.deaths.df1)
summary(model1) #0.3035
hist(model1$residuals, main = "Model residuals", col = "grey")
plot(model1$residuals ~ model1$fitted.values, xlab = "fitted values", ylab = "residuals")
qqnorm(model1$residuals, xlab = "Theoretical Quantiles: model1 residuals" )
qqline(model1$residuals, col = 2) ## red color

# relative importance of variables
library(relaimpo)
calc.relimp(model1, type = c("lmg"), rela = TRUE)

# calculate variance inflation factor
library(car)
vif(model1)
sqrt(vif(model1)) > 2  # if > 2 vif too high

# There's a huge multicollinearity & not really sure which variables to drop
# Kaiser-Meyer-Olkin statistics: if overall MSA > 0.6, proceed to factor analysis
library(psych)
KMO(cor(covid.deaths.df1))
# Overall MSA =  0.73, greater than 0.6 so let's proceed with PCA
# Determine Number of Factors to Extract
library(nFactors)
# get eigenvalues: eigen() uses a correlation matrix
ev <- eigen(cor(covid.deaths.df1))
ev$values
# plot a scree plot of eigenvalues
plot(ev$values, type="b", col="blue", xlab="variables")

# calculate cumulative proportion of eigenvalue and plot
ev.sum<-0
for(i in 1:length(ev$value)){
  ev.sum<-ev.sum+ev$value[i]
}
ev.list1<-1:length(ev$value)
for(i in 1:length(ev$value)){
  ev.list1[i]=ev$value[i]/ev.sum
}
ev.list2<-1:length(ev$value)
ev.list2[1]<-ev.list1[1]
for(i in 2:length(ev$value)){
  ev.list2[i]=ev.list2[i-1]+ev.list1[i]
}
plot (ev.list2, type="b", col="red", xlab="number of components", ylab ="cumulative proportion")
# 4 components

# Varimax Rotated Principal Components
# retaining 'nFactors' components
library(GPArotation)

# principal() uses a data frame or matrix of correlations
fit <- principal(covid.deaths.df1, nfactors=4, rotate="varimax")
fit
# create variables to represent the rotated components
fit$scores
cd.df1.fit.data <- data.frame(fit$scores)
# check new variables are uncorrelated
cor.matrix <-cor(cd.df1.fit.data, method = "pearson")
cor.df <- as.data.frame(cor.matrix)
round(cor.df, 2)

# model from new variables derived from factor analysis
model1b <- lm (Total_Deaths_1000 ~ cd.df1.fit.data$RC1 + cd.df1.fit.data$RC2 + cd.df1.fit.data$RC3 + cd.df1.fit.data$RC4)
summary(model1b) #0.1852
calc.relimp(model1b, type = c("lmg"), rela = TRUE)

# calculate variance inflation factor
vif(model1b)
sqrt(vif(model1b)) > 2
# This model doesn't have collinearity (but also doesn't fit that well)

# use a step-wise approach to search for a best model
library(RcmdrMisc)
# forward/backward step-wise selection
model1c <- stepwise(model1) # default direction = "backward/forward"
summary(model1c) #0.3027
model1d <- stepwise(model1, direction = "forward/backward")
summary(model1d) #0.3046
model1e <- stepwise(model1, direction = "forward")
summary(model1e) #0.3022
model1f <- stepwise(model1, direction = "backward")
summary(model1f) #0.3027
calc.relimp(model1, type = c("lmg"), rela = TRUE)

anova(model1c, model1d, test = "F")

#----------Section 7--------
# Partial Correlation

library(ppcor)
# Check correlation again but combining all normal & non-normal variables
pairs.panels(covid.deaths.df, method = "pearson", hist.col = "grey", col = "blue", main = "Pearson")
pairs.panels(covid.deaths.df, method = "spearman", hist.col = "grey", col = "blue", main = "Spearman")

cor.matrix3 <- cor(covid.deaths.df, use = "pairwise.complete.obs", method = "pearson")
corrplot(cor.matrix3, type = "upper", tl.col = "black", tl.srt = 45)
cor.matrix4 <- cor(covid.deaths.df, use = "pairwise.complete.obs", method = "spearman")
corrplot(cor.matrix4, type = "upper", tl.col = "black", tl.srt = 45)
# From corplots, it is evident that internal correlation among same themed variables is very high
# We'll choose variables that are either significant or correlated with dependent variable
age_theme <- c("X15_and_under", "X35_to_49", "X50_to_74", "X75_and_over")
health_theme <- c("Very_good_health", "Fair_health", "Bad_health", "Very_bad_health")
ethnicity_theme <- c("Asian_All", "Black_All", "Other")
deprivation_theme <- c("HH_no_depr", "HH_Depr_1_dim", "HH_Depr_2_dim", "HH_Depr_3_dim")
people_density_theme <- c("x_gt_1p5_to_2")

themes <- c("age_theme", "health_theme", "ethnicity_theme", "deprivation_theme", "people_density_theme")

pcor.res <- data.frame()
temp.res <- data.frame("var1" = character(), "controlled_vars" = character(), "coeff" = numeric(), "p-val" = numeric())
for (theme in themes) {
  for (var1 in get(theme)) {
    temp <- list()
    pcor.res <- pcor.test(Total_Deaths_1000, covid.deaths.df1[, var1], matrix(0, nrow = length(Total_Deaths_1000), ncol = 1))
    pcor.res.p <- ifelse(pcor.res$p.value < 1e-4,
                         format(pcor.res$p.value, scientific = TRUE),
                         format(pcor.res$p.value, scientific = FALSE))
    temp.res[nrow(temp.res)+1, ] <- list(var1, paste(unlist(temp), collapse = ", "), pcor.res$estimate, pcor.res.p)
    for (var2 in get(theme)) {
      if (var1 != var2) {
        pcor.res <- pcor.test(Total_Deaths_1000, covid.deaths.df1[, var1], covid.deaths.df1[, var2])
        pcor.res.p <- ifelse(pcor.res$p.value < 1e-4,
                             format(pcor.res$p.value, scientific = TRUE),
                             format(pcor.res$p.value, scientific = FALSE))
        temp.res[nrow(temp.res)+1, ] <- list(var1, var2, pcor.res$estimate, pcor.res.p)
        
        temp <- c(temp, var2)
        if (length(temp) == 1){
          next
        }
        pcor.res <- pcor.test(Total_Deaths_1000, covid.deaths.df1[, var1], covid.deaths.df1[, unlist(temp)])
        pcor.res.p <- ifelse(pcor.res$p.value < 1e-4,
                             format(pcor.res$p.value, scientific = TRUE),
                             format(pcor.res$p.value, scientific = FALSE))
        temp.res[nrow(temp.res)+1, ] <- list(var1, paste(unlist(temp), collapse = ", "), pcor.res$estimate, pcor.res.p)
      }
    }
  }
}

# Interesting things popped when checked for collinearity with sibling variables (as in variables under same parent theme)
# Pairs to be noticed & comment about them
# age_theme: X15_and_under shows improved correlation when all its sibling variables are controlled
#         It looked like either X15_and_under has to be removed or all other sibling variables have to be removed to get rid of partial correlation in that theme
#         (X15_and_under alone or variables except it) from age theme combined with all other variables explains relation with death better
# health_theme: Very good & Fair health can be removed as their effect on dependent variable while controlling Bad & Very bad health is less, and it is more when checked vice-versa
# ethnicity_theme: variables from these theme are chosen only because their significance, so checking correlation & doing partial correlation analysis is waste
# deprivation_theme: HH_no_depr & hh_2_dim_depr shows better correlation even when other sibling variables are controlled but vice-versa is not true
# people_density_theme: this can be ignored

covid.deaths.pcor1 <- covid.deaths[,c("X15_and_under",
                                     "Bad_health", "Very_bad_health",
                                     "Asian_All", "Black_All", "Other",
                                     "HH_no_depr", "HH_Depr_2_dim",
                                     "x_gt_1p5_to_2")]
covid.deaths.pcor2 <- covid.deaths[,c("X35_to_49", "X50_to_74", "X75_and_over",
                                      "Bad_health", "Very_bad_health",
                                      "Asian_All", "Black_All", "Other",
                                      "HH_no_depr", "HH_Depr_2_dim",
                                      "x_gt_1p5_to_2")]
model2a <- lm(Total_Deaths_1000 ~ ., data = covid.deaths.pcor1)
summary(model2a) #0.2021
model2b <- lm(Total_Deaths_1000 ~ ., data = covid.deaths.pcor2)
summary(model2b) #0.3079
hist(model2b$residuals, main = "Model residuals", col = "grey")
plot(model2b$residuals ~ model2b$fitted.values, xlab = "fitted values", ylab = "residuals")
qqnorm(model2b$residuals, xlab = "Theoretical Quantiles: model2b residuals" )
qqline(model2b$residuals, col = 2) ## red color
calc.relimp(model2b, type = c("lmg"), rela = TRUE)

# test whether model1 and model2b are significantly different using F test
anova(model2b, model1, test = "F")

# calculate variance inflation factor
vif(model2b)
sqrt(vif(model2b)) > 2  # if > 2 vif too high

# There's a huge multicollinearity & not really sure which variables to drop
# Kaiser-Meyer-Olkin statistics: if overall MSA > 0.6, proceed to factor analysis

KMO(cor(covid.deaths.pcor2))
# Overall MSA =  0.8, greater than 0.6 so let's proceed with PCA

ev <- eigen(cor(covid.deaths.pcor2))
# plot a scree plot of eigenvalues
plot(ev$values, type="b", col="blue", xlab="variables")

# calculate cumulative proportion of eigenvalue and plot
ev.sum<-0
for(i in 1:length(ev$value)){
  ev.sum<-ev.sum+ev$value[i]
}
ev.list1<-1:length(ev$value)
for(i in 1:length(ev$value)){
  ev.list1[i]=ev$value[i]/ev.sum
}
ev.list2<-1:length(ev$value)
ev.list2[1]<-ev.list1[1]
for(i in 2:length(ev$value)){
  ev.list2[i]=ev.list2[i-1]+ev.list1[i]
}
plot (ev.list2, type="b", col="red", xlab="number of components", ylab ="cumulative proportion")
# 3 components

fit <- principal(covid.deaths.pcor2, nfactors=3, rotate="varimax")
fit
# create variables to represent the rotated components
fit$scores
cd.pcor2.fit.data <- data.frame(fit$scores)
# check new variables are uncorrelated
cor.matrix <-cor(cd.pcor2.fit.data, method = "pearson")
cor.df <- as.data.frame(cor.matrix)
round(cor.df, 2)

# model from new variables derived from factor analysis
model2b_1 <- lm (Total_Deaths_1000 ~ cd.pcor2.fit.data$RC1 + cd.pcor2.fit.data$RC2 + cd.pcor2.fit.data$RC3)
summary(model2b_1) #0.1622

# calculate variance inflation factor
vif(model2b_1)
sqrt(vif(model2b_1)) > 2
# This model doesn't have collinearity (but also doesn't fit that well)

# use a stepwise approach to search for a best model
# forward/backward stepwise selection
model2b_2 <- stepwise(model2b)
summary(model2b_2) #0.291
model2b_3 <- stepwise(model2b, direction = "forward")
summary(model2b_3) #0.3004
#model2b_4 <- stepwise(model2b, direction = "forward/backward")
#summary(model2b_4) #0.291
#model2b_5 <- stepwise(model2b, direction = "backward")
#summary(model2b_5) #0.291

#----------Section 8--------
# Thematic PCA for dimensional reduction 

# PCA on Age 
age_df <- covid.deaths[,age_theme]
cor(age_df)
KMO(cor(age_df)) #0.77
# Learning from previous partial correlation analysis we'll remove variable "X15_and_under"
age_theme2 <- c("X35_to_49", "X50_to_74", "X75_and_over")
age_df2 <- covid.deaths[,age_theme2]
cor(age_df2)
KMO(cor(age_df2)) #0.72
ev <- eigen(cor(age_df2))
plot(ev$values, type="b", col="blue", xlab="variables")
# calculate cumulative proportion of eigenvalue and plot
ev.sum<-0
for(i in 1:length(ev$value)){
  ev.sum<-ev.sum+ev$value[i]
}
ev.list1<-1:length(ev$value)
for(i in 1:length(ev$value)){
  ev.list1[i]=ev$value[i]/ev.sum
}
ev.list2<-1:length(ev$value)
ev.list2[1]<-ev.list1[1]
for(i in 2:length(ev$value)){
  ev.list2[i]=ev.list2[i-1]+ev.list1[i]
}
plot (ev.list2, type="b", col="red", xlab="number of components", ylab ="cumulative proportion")
# 2 components
fit <- principal(age_df2, nfactors=2, rotate="varimax")
fit
# create variables to represent the rotated components
age_df2.fit.data <- data.frame(fit$scores)
colnames(age_df2.fit.data) <- c("Age_pca_1", "Age_pca_2")
# Age_pca_1 = People above age 50 (or) Senior People
# Age_pca_2 = 35 to 49
# check new variables are uncorrelated
cor(age_df2.fit.data, method = "pearson")

# PCA on General health 
health_df <- covid.deaths[,health_theme]
cor(health_df)
KMO(cor(health_df)) #0.69
ev <- eigen(cor(health_df))
plot(ev$values, type="b", col="blue", xlab="variables")
# calculate cumulative proportion of eigenvalue and plot
ev.sum<-0
for(i in 1:length(ev$value)){
  ev.sum<-ev.sum+ev$value[i]
}
ev.list1<-1:length(ev$value)
for(i in 1:length(ev$value)){
  ev.list1[i]=ev$value[i]/ev.sum
}
ev.list2<-1:length(ev$value)
ev.list2[1]<-ev.list1[1]
for(i in 2:length(ev$value)){
  ev.list2[i]=ev.list2[i-1]+ev.list1[i]
}
plot (ev.list2, type="b", col="red", xlab="number of components", ylab ="cumulative proportion")
# 2 components
fit <- principal(health_df, nfactors=2, rotate="varimax")
fit
# create variables to represent the rotated components
health_df.fit.data <- data.frame(fit$scores)
colnames(health_df.fit.data) <- c("Health_pca_1", "Health_pca_2")
# Health_pca_1 = Very_good_health & Fair_health
# Health_pca_2 = Bad_health & Very_bad_health
# check new variables are uncorrelated
cor(health_df.fit.data, method = "pearson")

#PCA on Ethnicity theme
ethnic_df <- covid.deaths[,ethnicity_theme]
cor(ethnic_df)
KMO(cor(ethnic_df)) #0.68
ev <- eigen(cor(ethnic_df))
plot(ev$values, type="b", col="blue", xlab="variables")
# calculate cumulative proportion of eigenvalue and plot
ev.sum<-0
for(i in 1:length(ev$value)){
  ev.sum<-ev.sum+ev$value[i]
}
ev.list1<-1:length(ev$value)
for(i in 1:length(ev$value)){
  ev.list1[i]=ev$value[i]/ev.sum
}
ev.list2<-1:length(ev$value)
ev.list2[1]<-ev.list1[1]
for(i in 2:length(ev$value)){
  ev.list2[i]=ev.list2[i-1]+ev.list1[i]
}
plot (ev.list2, type="b", col="red", xlab="number of components", ylab ="cumulative proportion")
# 2 components
fit <- principal(ethnic_df, nfactors=2, rotate="varimax")
fit
# create variables to represent the rotated components
ethnic_df.fit.data <- data.frame(fit$scores)
colnames(ethnic_df.fit.data) <- c("Ethnic_pca_1", "Ethnic_pca_2")
# Ethnic_pca_1 = Black_All & Other
# Ethnic_pca_2 = Asian_all
# check new variables are uncorrelated
cor(ethnic_df.fit.data, method = "pearson")

#PCA on HH_deprivation theme
depr_df <- covid.deaths[,deprivation_theme]
cor(depr_df)
KMO(cor(depr_df)) #0.4

covid.deaths.theme.pca <- data.frame(c(age_df2.fit.data, health_df.fit.data, ethnic_df.fit.data, depr_df, data.frame(x_gt_1p5_to_2, Total_Deaths_1000)))

pairs.panels(covid.deaths.theme.pca, method = "pearson", hist.col = "grey", col = "blue", main = "Pearson")

# Exclude dependent variable 
covid.deaths.theme.pca.df <- covid.deaths.theme.pca[, setdiff(names(covid.deaths.theme.pca), "Total_Deaths_1000")]
model3 <- lm(Total_Deaths_1000 ~ ., data = covid.deaths.theme.pca.df)
summary(model3) #0.309
calc.relimp(model3, type = c("lmg"), rela = TRUE)

# calculate variance inflation factor
vif(model3)
sqrt(vif(model3)) > 2

KMO(cor(covid.deaths.theme.pca.df))

# MSA = 0.52 , below 0.6 so can't do cumulative factor analysis

# let's see if we can remove variables using partial correlation
temp2.res <- data.frame("var1" = character(), "controlling_var" = character(), "coeff" = numeric(), "p-val" = numeric())
names_list <- as.list(names(covid.deaths.theme.pca.df))
for (var1 in names_list) {
  names_list[[1]] <- NULL
  for (var2 in names_list) {
    pcor.res <- pcor.test(Total_Deaths_1000, covid.deaths.theme.pca.df[, var1], covid.deaths.theme.pca.df[,var2])
    pcor.res.p <- ifelse(pcor.res$p.value < 5e-2,
                         format(pcor.res$p.value, scientific = TRUE),
                         format(pcor.res$p.value, scientific = FALSE))
    temp2.res[nrow(temp2.res)+1, ] <- list(var1, var2, pcor.res$estimate, pcor.res.p)
    pcor.res <- pcor.test(Total_Deaths_1000, covid.deaths.theme.pca.df[, var2], covid.deaths.theme.pca.df[,var1])
    pcor.res.p <- ifelse(pcor.res$p.value < 5e-2,
                         format(pcor.res$p.value, scientific = TRUE),
                         format(pcor.res$p.value, scientific = FALSE))
    temp2.res[nrow(temp2.res)+1, ] <- list(var2, var1, pcor.res$estimate, pcor.res.p)
  }
}
# we can remove Age_pca_2, HH_Depr_1_dim 
covid.deaths.theme.pca.df2 <- covid.deaths.theme.pca.df[, setdiff(names(covid.deaths.theme.pca.df), c("Age_pca_2", "HH_Depr_1_dim"))]
model3b <- lm(Total_Deaths_1000 ~ ., data = covid.deaths.theme.pca.df2)
summary(model3b) #0.3128
calc.relimp(model3b, type = c("lmg"), rela = TRUE)

anova(model3, model3b, test = "F")

# calculate variance inflation factor
vif(model3b)
sqrt(vif(model3b)) > 2

KMO(cor(covid.deaths.theme.pca.df2)) #0.51

# use a stepwise approach to search for a best model
# forward/backward stepwise selection
model3c <- stepwise(model3)
summary(model3c) #0.2925
model3d <- stepwise(model3, direction = "forward")
summary(model3d) #0.2834
#model3e <- stepwise(model3, direction = "forward/backward")
#summary(model3e) #0.2834
#model3f <- stepwise(model3, direction = "backward")
#summary(model3f) #0.2925
model3g <- stepwise(model3b)
summary(model3g) #0.2834
model3f <- stepwise(model3b, direction = "forward")
summary(model3f) #0.2834

#----------Section 9--------
# Thematic-PCA on raw data

# Check KMO on themes
all_age_theme <- c("X15_and_under", "X16_to_24", "X25_to_34", "X35_to_49", "X50_to_74", "X75_and_over")
all_health_theme <- c("Very_good_health", "Good_health", "Fair_health", "Bad_health", "Very_bad_health")
all_ethnicity_theme <- c("Asian_All", "Black_All", "Mixed", "White", "Other")
all_deprivation_theme <- c("HH_no_depr", "HH_Depr_1_dim", "HH_Depr_2_dim", "HH_Depr_3_dim", "HH_Depr_4_dim")
all_people_density_theme <- c("x_lt_0p5", "x_0p5_to_lt_1", "x_eq_1", "x_gt_1_to_1p5", "x_gt_1p5_to_2", "x_gt_2")

all_age_df <- covid.deaths[,all_age_theme]
KMO(all_age_df) #0.28
all_health_df <- covid.deaths[,all_health_theme]
KMO(all_health_df) #0.35
all_ethnicity_df <- covid.deaths[,all_ethnicity_theme]
KMO(all_ethnicity_df) #0.36
all_deprivation_df <- covid.deaths[,all_deprivation_theme]
KMO(all_deprivation_df) #0.34
all_people_density_df <- covid.deaths[,all_people_density_theme]
KMO(all_people_density_df) #0.67

ev <- eigen(cor(all_people_density_df))
plot(ev$values, type="b", col="blue", xlab="variables")
# calculate cumulative proportion of eigenvalue and plot
ev.sum<-0
for(i in 1:length(ev$value)){
  ev.sum<-ev.sum+ev$value[i]
}
ev.list1<-1:length(ev$value)
for(i in 1:length(ev$value)){
  ev.list1[i]=ev$value[i]/ev.sum
}
ev.list2<-1:length(ev$value)
ev.list2[1]<-ev.list1[1]
for(i in 2:length(ev$value)){
  ev.list2[i]=ev.list2[i-1]+ev.list1[i]
}
plot (ev.list2, type="b", col="red", xlab="number of components", ylab ="cumulative proportion")
# 2 components
fit <- principal(all_people_density_df, nfactors=2, rotate="varimax")
fit
all_people_density_df.fit.data <- data.frame(fit$scores)
colnames(all_people_density_df.fit.data) <- c("people_density_pca_1", "people_density_pca_2")
# check new variables are uncorrelated
cor(all_people_density_df.fit.data, method = "pearson")

covid.deaths.df2 <- covid.deaths[,setdiff(names(covid.deaths), c(all_people_density_theme, "Total_Deaths_1000"))]
covid.deaths.df2 <- cbind(covid.deaths.df2, all_people_density_df.fit.data)

model4 <- lm(Total_Deaths_1000 ~ ., data = covid.deaths.df2)
summary(model4) #0.3045
calc.relimp(model4, type = c("lmg"), rela = TRUE)

anova(model3b, model4, test = "F")
# calculate variance inflation factor
vif(model4)
sqrt(vif(model4)) > 2  # if > 2 vif too high

# There's a huge multicollinearity & not really sure which variables to drop
# Kaiser-Meyer-Olkin statistics: if overall MSA > 0.6, proceed to factor analysis

KMO(cor(covid.deaths.df2))

# Overall MSA =  0.58, close to 0.6 so let's see if PCA helps
ev <- eigen(cor(covid.deaths.df2))
plot(ev$values, type="b", col="blue", xlab="variables")

# calculate cumulative proportion of eigenvalue and plot
ev.sum<-0
for(i in 1:length(ev$value)){
  ev.sum<-ev.sum+ev$value[i]
}
ev.list1<-1:length(ev$value)
for(i in 1:length(ev$value)){
  ev.list1[i]=ev$value[i]/ev.sum
}
ev.list2<-1:length(ev$value)
ev.list2[1]<-ev.list1[1]
for(i in 2:length(ev$value)){
  ev.list2[i]=ev.list2[i-1]+ev.list1[i]
}
plot (ev.list2, type="b", col="red", xlab="number of components", ylab ="cumulative proportion")
# 5 components

# Varimax Rotated Principal Components
# retaining 'nFactors' components
fit <- principal(covid.deaths.df2, nfactors=5, rotate="varimax")
fit

# create variables to represent the rotated components
cd.df2.fit.data <- data.frame(fit$scores)
# check new variables are uncorrelated
cor(cd.df2.fit.data, method = "pearson")

model4b <- lm (Total_Deaths_1000 ~ cd.df2.fit.data$RC1 + cd.df2.fit.data$RC2 + cd.df2.fit.data$RC3 + cd.df2.fit.data$RC4 + cd.df2.fit.data$RC5)
summary(model4b) #0.2536

# Examining loading data in factor analysis of covid.deaths.df2
# Variables X35_to_49, X75_and_over have complexity > 3
myvars <- names(covid.deaths.df2) %in% c("X35_to_49", "X75_and_over")
covid.deaths.df3 <- covid.deaths.df2[!myvars]
rm(myvars)

KMO(cor(covid.deaths.df3))

# Overall MSA =  0.59, close to 0.6 so let's proceed with PCA
ev <- eigen(cor(covid.deaths.df3))
plot(ev$values, type="b", col="blue", xlab="variables")
# calculate cumulative proportion of eigenvalue and plot
ev.sum<-0
for(i in 1:length(ev$value)){
  ev.sum<-ev.sum+ev$value[i]
}
ev.list1<-1:length(ev$value)
for(i in 1:length(ev$value)){
  ev.list1[i]=ev$value[i]/ev.sum
}
ev.list2<-1:length(ev$value)
ev.list2[1]<-ev.list1[1]
for(i in 2:length(ev$value)){
  ev.list2[i]=ev.list2[i-1]+ev.list1[i]
}
plot (ev.list2, type="b", col="red", xlab="number of components", ylab ="cumulative proportion")
# 5 components

fit <- principal(covid.deaths.df3, nfactors=5, rotate="varimax")
fit

# create variables to represent the rotated components
cd.df3.fit.data <- data.frame(fit$scores)
# check new variables are uncorrelated
cor(cd.df3.fit.data, method = "pearson")

# model from new variables derived from factor analysis
model4c <- lm (Total_Deaths_1000 ~ cd.df3.fit.data$RC1 + cd.df3.fit.data$RC2 + cd.df3.fit.data$RC3 + cd.df3.fit.data$RC4 + cd.df3.fit.data$RC5)
summary(model4c) #0.2542

# calculate variance inflation factor
vif(model4c)
sqrt(vif(model4c)) > 2
# This model doesn't have collinearity (but also doesn't fit that well)

# use a stepwise approach to search for a best model
# forward/backward stepwise selection
model4d <- stepwise(model4)
summary(model4d) #0.3032
model4g <- stepwise(model4, direction = "forward")
summary(model4g) #0.2991
#model4h <- stepwise(model4, direction = "forward/backward")
#summary(model4h) #0.2991
#model4i <- stepwise(model4, direction = "backward")
#summary(model4i) #0.3032
model4e <- stepwise(model4b)
summary(model4e) #0.2559
model4f <- stepwise(model4c)
summary(model4f) #0.2557


#----------Section 10--------
# Regression on raw data

covid.deaths.all <- covid.deaths[,setdiff(names(covid.deaths), "Total_Deaths_1000")]
model5 <- lm(Total_Deaths_1000 ~ ., data = covid.deaths.all)
summary(model5) #0.3179
calc.relimp(model5, type = c("lmg"), rela = TRUE)

anova(model5, model1, test = "F")
anova(model5, model2b, test = "F")
anova(model5, model3, test = "F")
anova(model5, model4, test = "F")

vif(model5)
sqrt(vif(model5)) > 2  # if > 2 vif too high

# There's a huge multicollinearity & not really sure which variables to drop
# Kaiser-Meyer-Olkin statistics: if overall MSA > 0.6, proceed to factor analysis

KMO(cor(covid.deaths.all))

# Overall MSA =  0.66, greater than 0.6 so let's proceed with PCA
ev <- eigen(cor(covid.deaths.all))
ev$values
# plot a scree plot of eigenvalues
plot(ev$values, type="b", col="blue", xlab="variables")

# calculate cumulative proportion of eigenvalue and plot
ev.sum<-0
for(i in 1:length(ev$value)){
  ev.sum<-ev.sum+ev$value[i]
}
ev.list1<-1:length(ev$value)
for(i in 1:length(ev$value)){
  ev.list1[i]=ev$value[i]/ev.sum
}
ev.list2<-1:length(ev$value)
ev.list2[1]<-ev.list1[1]
for(i in 2:length(ev$value)){
  ev.list2[i]=ev.list2[i-1]+ev.list1[i]
}
plot (ev.list2, type="b", col="red", xlab="number of components", ylab ="cumulative proportion")
# 5 components

fit <- principal(covid.deaths.all, nfactors=5, rotate="varimax")
fit

# Examining loading data in factor analysis of covid.deaths.all
# Variables X35_to_49, X75_and_over have complexity > 3
myvars <- names(covid.deaths.all) %in% c("X35_to_49", "X75_and_over")
covid.deaths.all <- covid.deaths.all[!myvars]
rm(myvars)

KMO(cor(covid.deaths.all))

# Overall MSA =  0.67, more than 0.6 so let's proceed with PCA
ev <- eigen(cor(covid.deaths.all))
plot(ev$values, type="b", col="blue", xlab="variables")
# calculate cumulative proportion of eigenvalue and plot
ev.sum<-0
for(i in 1:length(ev$value)){
  ev.sum<-ev.sum+ev$value[i]
}
ev.list1<-1:length(ev$value)
for(i in 1:length(ev$value)){
  ev.list1[i]=ev$value[i]/ev.sum
}
ev.list2<-1:length(ev$value)
ev.list2[1]<-ev.list1[1]
for(i in 2:length(ev$value)){
  ev.list2[i]=ev.list2[i-1]+ev.list1[i]
}
plot (ev.list2, type="b", col="red", xlab="number of components", ylab ="cumulative proportion")
# 5 components

fit <- principal(covid.deaths.all, nfactors=5, rotate="varimax")
fit

# create variables to represent the rotated components
cd.all.fit.data <- data.frame(fit$scores)
# check new variables are uncorrelated
cor(cd.all.fit.data, method = "pearson")

# model from new variables derived from factor analysis
model5b <- lm (Total_Deaths_1000 ~ cd.all.fit.data$RC1 + cd.all.fit.data$RC2 + cd.all.fit.data$RC3 + cd.all.fit.data$RC4 + cd.all.fit.data$RC5)
summary(model5b) #0.2524

# use a stepwise approach to search for a best model
# forward/backward stepwise selection
model5c <- stepwise(model5)
summary(model5c) #0.3244
model5e <- stepwise(model5, direction = "forward")
summary(model5e) #0.3148
#model5f <- stepwise(model5, direction = "forward/backward")
#summary(model5f) #0.3148
#model5g <- stepwise(model5, direction = "backward")
#summary(model5g) #0.3244
model5d <- stepwise(model5b)
summary(model5d) #0.2537

anova(model5c, model1c, test = "F")

#----------Section 11----------------

detach(covid.deaths)
# remove all variables from the environment
rm(list=ls())
