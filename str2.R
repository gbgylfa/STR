library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(reshape2)
library(memisc)
library(lme4)
library(texreg)
library(afex)

data <- read.csv('/Users/Duna/Documents/Ling/STR/data.csv')
head(data)

#### just initials######

init <- read.csv('/Users/Duna/Documents/Ling/STR/toplot.csv')
init_str <- subset(init,init$segment == "STR")
init_IHELP_str <- subset(init_str,init_str$Corpus == "IHELP")



plot_init_normalized_by_dob <- ggplot(data = init, aes(x = DOB, y = normalizedmean, colour = segment)) +
  geom_point(size=5)
print(plot_init_normalized_by_dob)

plot_init_by_dob <- ggplot(data = init, aes(x = DOB, y = meancog, colour = segment)) +
  geom_point(size=5)
print(plot_init_by_dob)

plot_init_by_dob_str <- ggplot(data = init_str, aes(x = DOB, y = normalizedmean, colour = segment)) +
  geom_point(size=5)
print(plot_init_by_dob_str)

plot_init_by_dob_IHELP_str <- ggplot(data = init_IHELP_str, aes(x = DOB, y = normalizedmean, colour = HS_code)) +
  geom_point(size=5)
print(plot_init_by_dob_IHELP_str)

plot_init_by_bim_IHELP_str <- ggplot(data = init_IHELP_str, aes(x = Phila_bimodality, y = normalizedmean, colour = HS_code)) +
  geom_point(size=5)
print(plot_init_by_bim_IHELP_str)

plot_init_by_nas_IHELP_str <- ggplot(data = init_IHELP_str, aes(x = nasal_bimodality, y = normalizedmean, colour = HS_code)) +
  geom_point(size=5)
print(plot_init_by_nas_IHELP_str)

plot_init_by_dob_str_income <- ggplot(data = init_str, aes(x = DOB, y = normalizedmean, colour = Income)) +
  geom_point(size=5)
print(plot_init_by_dob_str_income)

plot_init_by_income_str <- ggplot(data = init_str, aes(x = Income, y = normalizedmean, colour = DOB)) +
  geom_point(size=5)
print(plot_init_by_income_str)

plot_init_by_dob_str_ed <- ggplot(data = init_str, aes(x = DOB, y = normalizedmean, colour = EdYears)) +
  geom_point(size=5)
print(plot_init_by_dob_str_ed)

plot_init_by_ed_str <- ggplot(data = init_str, aes(x = EdYears, y = normalizedmean, colour = DOB)) +
  geom_point(size=5)
print(plot_init_by_ed_str)

plot_init_by_dob_str_sex <- ggplot(data = init_str, aes(x = DOB, y = normalizedmean, colour = Sex)) +
  geom_point(size=5)
print(plot_init_by_dob_str_sex)



######now all that for all of it ##########

all <- read.csv('/Users/Duna/Documents/Ling/STR/toplot_all.csv')
diff <-  read.csv('/Users/Duna/Documents/Ling/STR/toplot_diff.csv')
all_str <- subset(all,segment=="STR")
all_IHELP_str <- subset(all_str,Corpus=="IHELP")


plot_all_normalized_by_dob <- ggplot(data = all, aes(x = DOB, y = normalizedmean, colour = segment)) +
  geom_point(size=5) 
print(plot_all_normalized_by_dob)

# plot_all_normalized_by_dob_smooth <- ggplot(data = all, aes(x = DOB, y = normalizedmean, colour = segment)) +
#   geom_point(size=5) +
#   geom_smooth(data = all[all$segment=="STR",], method = 'auto') +
#   geom_smooth(data = all[all$segment=="SH",], method = 'auto') +
#   geom_smooth(data = all[all$segment=="S",], method = 'auto')
# print(plot_all_normalized_by_dob_smooth)

plot_all_normalized_by_dob_not_s <- ggplot(data = all[all$segment != "S",], aes(x = DOB, y = normalizedmean, colour = segment)) +
  geom_point(size=5) +
  scale_colour_manual(values=c("green", "blue"))
print(plot_all_normalized_by_dob_not_s)

plot_all_by_dob <- ggplot(data = all, aes(x = DOB, y = meancog, colour = segment)) +
  geom_point(size=5)
print(plot_all_by_dob)

plot_all_by_dob_str <- ggplot(data = all_str, aes(x = DOB, y = normalizedmean, colour = segment)) +
  geom_point(size=5)
print(plot_all_by_dob_str)

plot_all_by_dob_IHELP_str <- ggplot(data = all_IHELP_str, aes(x = DOB, y = normalizedmean, colour = HS_code)) +
  geom_point(size=5)
print(plot_all_by_dob_IHELP_str)

plot_all_by_bim_IHELP_str <- ggplot(data = all_IHELP_str, aes(x = Phila_bimodality, y = normalizedmean, colour = HS_code)) +
  geom_point(size=5)
print(plot_all_by_bim_IHELP_str)

plot_all_by_nas_IHELP_str <- ggplot(data = all_IHELP_str, aes(x = nasal_bimodality, y = normalizedmean, colour = HS_code)) +
  geom_point(size=5)
print(plot_all_by_nas_IHELP_str)

plot_all_by_dob_str_income <- ggplot(data = all_str, aes(x = DOB, y = normalizedmean, colour = Income)) +
  geom_point(size=5)
print(plot_all_by_dob_str_income)

plot_all_by_income_str <- ggplot(data = all_str, aes(x = Income, y = normalizedmean, colour = DOB)) +
  geom_point(size=5)
print(plot_all_by_income_str)

plot_all_by_dob_str_ed <- ggplot(data = all_str, aes(x = DOB, y = normalizedmean, colour = EdYears)) +
  geom_point(size=5)
print(plot_all_by_dob_str_ed)

plot_all_by_ed_str <- ggplot(data = all_str, aes(x = EdYears, y = normalizedmean, colour = DOB)) +
  geom_point(size=5)
print(plot_all_by_ed_str)

plot_all_by_dob_str_sex <- ggplot(data = all_str, aes(x = DOB, y = normalizedmean, colour = Sex)) +
  geom_point(size=5)
print(plot_all_by_dob_str_sex)

plot_all_by_dob_diff <- ggplot(data = diff, aes(x = DOB, y = diff)) +
  geom_point(size=5)
print(plot_all_by_dob_diff)



###### stats #########
data$following_frontness <- relevel(data$following_frontness, ref="non-front")
data$following_height <- relevel(data$following_height, ref="non-high")
data$following_hf <- relevel(data$following_hf, ref="not")

data_str <- subset(data,data$segment == "STR")

data_str$DOB <- scale(data_str$DOB, center = T, scale = T)

word_freq <- read.csv('/Users/Duna/Documents/Ling/STR/word_freq.csv')
word_freq$word <- toupper(word_freq$word)
word_freq$freq <- as.numeric(word_freq$freqcount)
word_freq$freq <- scale(word_freq$freq, center = T, scale = T)

data_str = merge(data_str,subset(word_freq, freqcount != "#N/A"), by = "word")


data_fem<- subset(data_str,data_str$Sex == "f")
data_IHELP <- subset(data_str, data_str$Corpus == "IHELP")
data_PNC <- subset(data_str, data_str$Corpus == "PNC")
data_income <- subset(data_str, data_str$Income != 'NA')
data_income$Income = scale(data_income$Income,center = T, scale = T)
data_ed <- subset(data_str, data_str$EdYears != 'NA')
data_bim <- subset(data_IHELP, data_IHELP$nasal_bimodality != "NA")
                   
mod1<-lmer(data_str$normalizedmean ~   data_str$DOB * data_str$position + data_str$Sex  + data_str$following_hf + data_str$street  + (1|data_str$subject), REML = F)
summary(mod1)

mod1.1<-lmer(data_str$normalizedmean ~ data_str$DOB * data_str$Sex + data_str$street+ (1|data_str$subject), REML = F)
summary(mod1.1)

mod1.2<-mixed(normalizedmean ~   DOB + position + Sex  + following_hf + freq +street  + (1|subject), data = data_str)
summary(mod1.2)


mod2<-lmer(data_str$normalizedmean ~ data_str$DOB * data_str$Sex + data_str$following +  data_str$word + (1|data_str$subject), REML = F)
summary(mod2)

mod3<-lmer(data_fem$normalizedmean ~ data_fem$DOB  +  data_fem$street + data_fem$following +  (1|data_fem$subject), REML = F)
summary(mod3)

mod4<-lmer(data_IHELP$normalizedmean ~ data_IHELP$DOB * data_IHELP$Sex + data_IHELP$following +  data_IHELP$street + (1|data_IHELP$subject), REML = F)
summary(mod4)


mod5<-lmer(data_PNC$normalizedmean ~ data_PNC$DOB * data_PNC$Sex + data_PNC$following +  data_PNC$street + (1|data_PNC$subject), REML = F)
summary(mod5)

mod6<-lmer(data_income$normalizedmean ~ data_income$Income + data_income$DOB + data_income$Sex + data_income$following +  data_income$street + (1|data_income$subject), REML = F)
summary(mod6)


mod8<-lmer(data_ed$normalizedmean ~ data_ed$EdYears + data_ed$DOB * data_ed$Sex + data_ed$following +  data_ed$street + (1|data_ed$subject), REML = F)
summary(mod8)

mod9<-lmer(data_bim$normalizedmean ~ data_bim$Phila_bimodality + data_bim$DOB * data_bim$Sex + data_bim$following +  data_bim$street + (1|data_bim$subject), REML = F)
summary(mod9)

mod6.2<-lmer(normalizedmean ~   DOB + position + +Income + Sex  + following_hf + freq +street  + (1|subject), data = data_PNC)
summary(mod6.2)

modfull <- lmer(normalizedmean ~ DOB +  Sex + position  + following_frontness +
                  DOB*street  + freq + 
                  following_hf +
                  (1|subject),data = data_str, REML = F)
summary(modfull)

mixed(normalizedmean ~ DOB +  Sex +position  + following_frontness +
                  DOB*street  + freq + 
                  following_hf +
                  (1|subject),data = data_str, REML = F)

modwovowel <- lmer(normalizedmean ~   DOB *Sex + DOB*position  + DOB *street  + DOB *freq + 
                  (1|subject),data = data_str, REML = F)

modwostreet <- lmer(normalizedmean ~   DOB *Sex + DOB*position  + 
                  DOB * following_frontness +  DOB *freq + 
                  following_hf +
                  (1|subject),data = data_str, REML = F)

modwopos <- lmer(normalizedmean ~   DOB *Sex + 
                  DOB * following_frontness + DOB *street  + DOB *freq + 
                  following_hf +
                  (1|subject),data = data_str, REML = F)

modwopossex <-lmer(normalizedmean ~  
                  DOB * following_frontness + DOB *street  + DOB *freq + 
                  following_hf +
                  (1|subject),data = data_str, REML = F)

modwoposhf <-lmer(normalizedmean ~   DOB *Sex +   
                  DOB * following_frontness + DOB *street  + DOB *freq +
                  (1|subject),data = data_str, REML = F)

modwoposfront <-lmer(normalizedmean ~  DOB *Sex +    DOB *street  + DOB *freq + 
                  following_hf +
                  (1|subject),data = data_str, REML = F)

modwoposfreq <- lmer(normalizedmean ~   DOB *Sex + 
                  DOB * following_frontness +DOB * street   + 
                  following_hf +
                  (1|subject),data = data_str, REML = F)

modwoposfreqmodwoposfreqstint <- lmer(normalizedmean ~   DOB *Sex + 
                  DOB * following_frontness + street   + 
                  following_hf +
                  (1|subject),data = data_str, REML = F)

modwoposfreqmodwoposfreqvowint <- lmer(normalizedmean ~  DOB * Sex + 
                   following_frontness + DOB *street   + 
                  following_hf +
                  (1|subject),data = data_str, REML = F)

#best?
modwoposfreqmodwoposfreqvowintsexint <- lmer(normalizedmean ~   Sex + 
                   following_frontness + DOB *street   + 
                  following_hf +
                  (1|subject),data = data_str, REML = F)

modwoposfreqmodwoposhffreqvowintsexint <- lmer(normalizedmean ~   Sex + 
                   following_frontness + DOB *street +
                  (1|subject),data = data_str, REML = F)


rsquared.glmm(list(mod1,modfull,modwovowel,modwostreet,modwopos,modwopossex,modwoposhf,modwoposfront,modwoposfreqmodwoposfreqstint,
                   modwoposfreqmodwoposfreqvowint,modwoposfreqmodwoposhffreqvowintsexint))

modfinal <- lmer(normalizedmean ~   DOB *street   + 
                  following_hf +
                  (1|subject), data = data_str)
summary(modfinal)

mixed(normalizedmean ~   DOB *street   + 
                  following_hf +
                  (1|subject), data = data_str)

###### print ####
texreg(modfinal)
texreg(modfinal, t= TRUE)




#####income stuff#####

modfull <- lmer(normalizedmean ~   DOB + Income +Sex +position  + 
                  following_frontness + street  + freq + 
                  following_hf +
                  (1|subject),data = data_income, REML = F)

summary(modfull)


mixed(normalizedmean ~  DOB + Income +Sex +position  + 
                  following_frontness + street  + freq + 
                  following_hf +
                  (1|subject),data = data_income, REML = F)

summary(modfull)

modwovowel <- lmer(normalizedmean ~ Income +  DOB *Sex + DOB*position  + DOB *street  + DOB *freq + 
                  (1|subject),data = data_income, REML = F)

modwostreetint <- lmer(normalizedmean ~ Income +  DOB *Sex + street +DOB*position  + 
                  DOB * following_frontness +  DOB *freq + 
                  following_hf +
                  (1|subject),data = data_income, REML = F)

modwoposint <- lmer(normalizedmean ~ Income +  DOB *Sex + street +position  + 
                  DOB * following_frontness +  DOB *freq + 
                  following_hf +
                  (1|subject),data = data_income, REML = F)

modwoposinthf <- lmer(normalizedmean ~ Income +  DOB *Sex + street +position  + 
                  DOB * following_hf+  DOB *freq + 
                  following_hf +
                  (1|subject),data = data_income, REML = F)

modwohfint <- lmer(normalizedmean ~ Income +  Sex + street +position  + 
                DOB + following_hf+  freq + 
                  (1|subject),data = data_income, REML = F)

summary(modwohfint)


modwofreqint <- lmer(normalizedmean ~ Income +  DOB *Sex + street +position  + 
                following_hf+  freq +
                  following_hf +
                  (1|subject),data = data_income, REML = F)


rsquared.glmm(list(modfull,modwovowel,modwostreetint, modwoposint,modwoposinthf,modwohfint,modwofreqint))



#rsquared function

#' R-squared and pseudo-rsquared for a list of (generalized) linear (mixed) models
#'
#' This function calls the generic \code{\link{r.squared}} function for each of the
#' models in the list and rbinds the outputs into one data frame
#'
#' @param a list of fitted (generalized) linear (mixed) model objects
#' @return a dataframe with one row per model, and "Class",
#'         "Family", "Marginal", "Conditional" and "AIC" columns
rsquared.glmm <- function(modlist) {
  # Iterate over each model in the list
  do.call(rbind, lapply(modlist, r.squared))
}
 
#' R-squared and pseudo-rsquared for (generalized) linear (mixed) models
#'
#' This generic function calculates the r squared and pseudo r-squared for
#' a variety of(generalized) linear (mixed) model fits.
#' Currently implemented for \code{\link{lm}}, \code{\link{lmerTest::merMod}},
#' and \code{\link{nlme::lme}} objects.
#' Implementing methods usually call \code{\link{.rsquared.glmm}}
#'
#' @param mdl a fitted (generalized) linear (mixed) model object
#' @return Implementing methods usually return a dataframe with "Class",
#'         "Family", "Marginal", "Conditional", and "AIC" columns
r.squared <- function(mdl){
  UseMethod("r.squared")
}
 
#' Marginal r-squared for lm objects
#'
#' This method uses r.squared from \code{\link{summary}} as the marginal.
#' Contrary to other \code{\link{r.squared}} methods, 
#' this one doesn't call \code{\link{.rsquared.glmm}}
#'
#' @param mdl an lm object (usually fit using \code{\link{lm}},
#' @return a dataframe with with "Class" = "lm", "Family" = "gaussian",
#'        "Marginal" = unadjusted r-squared, "Conditional" = NA, and "AIC" columns
r.squared.lm <- function(mdl){
  data.frame(Class=class(mdl), Family="gaussian", Link="identity",
             Marginal=summary(mdl)$r.squared,
             Conditional=NA, AIC=AIC(mdl))
}
 
#' Marginal and conditional r-squared for merMod objects
#'
#' This method extracts the variance for fixed and random effects, residuals,
#' and the fixed effects for the null model (in the case of Poisson family),
#' and calls \code{\link{.rsquared.glmm}}
#'
#' @param mdl an merMod model (usually fit using \code{\link{lme4::lmer}},
#'        \code{\link{lme4::glmer}}, \code{\link{lmerTest::lmer}},
#'        \code{\link{blme::blmer}}, \code{\link{blme::bglmer}}, etc)
r.squared.merMod <- function(mdl){
  # Get variance of fixed effects by multiplying coefficients by design matrix
  VarF <- var(as.vector(lme4::fixef(mdl) %*% t(mdl@pp$X)))
  # Get variance of random effects by extracting variance components
  # Omit random effects at the observation level, variance is factored in later
  VarRand <- sum(
    sapply(
      VarCorr(mdl)[!sapply(unique(unlist(strsplit(names(ranef(mdl)),":|/"))), function(l) length(unique(mdl@frame[,l])) == nrow(mdl@frame))],
      function(Sigma) {
        X <- model.matrix(mdl)
        Z <- X[,rownames(Sigma)]
        sum(diag(Z %*% Sigma %*% t(Z)))/nrow(X) } ) )
  # Get the dispersion variance
  VarDisp <- unlist(VarCorr(mdl)[sapply(unique(unlist(strsplit(names(ranef(mdl)),":|/"))), function(l) length(unique(mdl@frame[,l])) == nrow(mdl@frame))])
  if(is.null(VarDisp)) VarDisp = 0 else VarDisp = VarDisp
  if(inherits(mdl, "lmerMod")){
    # Get residual variance
    VarResid <- attr(lme4::VarCorr(mdl), "sc")^2
    # Get ML model AIC
    mdl.aic <- AIC(update(mdl, REML=F))
    # Model family for lmer is gaussian
    family <- "gaussian"
    # Model link for lmer is identity
    link <- "identity"
  }
  else if(inherits(mdl, "glmerMod")){
    # Get the model summary
    mdl.summ <- summary(mdl)
    # Get the model's family, link and AIC
    family <- mdl.summ$family
    link <- mdl.summ$link
    mdl.aic <- AIC(mdl)
    # Pseudo-r-squared for poisson also requires the fixed effects of the null model
    if(family=="poisson") {
      # Get random effects names to generate null model
      rand.formula <- reformulate(sapply(findbars(formula(mdl)),
                                         function(x) paste0("(", deparse(x), ")")),
                                  response=".")
      # Generate null model (intercept and random effects only, no fixed effects)
      null.mdl <- update(mdl, rand.formula)
      # Get the fixed effects of the null model
      null.fixef <- as.numeric(lme4::fixef(null.mdl))
    }
  }
  # Call the internal function to do the pseudo r-squared calculations
  .rsquared.glmm(VarF, VarRand, VarResid, VarDisp, family = family, link = link,
                 mdl.aic = mdl.aic,
                 mdl.class = class(mdl),
                 null.fixef = null.fixef)
}
 
#' Marginal and conditional r-squared for lme objects
#'
#' This method extracts the variance for fixed and random effects,
#' as well as residuals, and calls \code{\link{.rsquared.glmm}}
#'
#' @param mdl an lme model (usually fit using \code{\link{nlme::lme}})
r.squared.lme <- function(mdl){
  # Get design matrix of fixed effects from model
  Fmat <- model.matrix(eval(mdl$call$fixed)[-2], mdl$data)
  # Get variance of fixed effects by multiplying coefficients by design matrix
  VarF <- var(as.vector(nlme::fixef(mdl) %*% t(Fmat)))
  # Get variance of random effects by extracting variance components
  VarRand <- sum(suppressWarnings(as.numeric(nlme::VarCorr(mdl)
                                             [rownames(nlme::VarCorr(mdl)) != "Residual",
                                              1])), na.rm=T)
  # Get residual variance
  VarResid <- as.numeric(nlme::VarCorr(mdl)[rownames(nlme::VarCorr(mdl))=="Residual", 1])
  # Call the internal function to do the pseudo r-squared calculations
  .rsquared.glmm(VarF, VarRand, VarResid, family = "gaussian", link = "identity",
                 mdl.aic = AIC(update(mdl, method="ML")),
                 mdl.class = class(mdl))
}
 
#' Marginal and conditional r-squared for glmm given fixed and random variances
#'
#' This function is based on Nakagawa and Schielzeth (2013). It returns the marginal
#' and conditional r-squared, as well as the AIC for each glmm.
#' Users should call the higher-level generic "r.squared", or implement a method for the
#' corresponding class to get varF, varRand and the family from the specific object
#'
#' @param varF Variance of fixed effects
#' @param varRand Variance of random effects
#' @param varResid Residual variance. Only necessary for "gaussian" family
#' @param family family of the glmm (currently works with gaussian, binomial and poisson)
#' @param link model link function. Working links are: gaussian: "identity" (default);
#'        binomial: "logit" (default), "probit"; poisson: "log" (default), "sqrt"
#' @param mdl.aic The model's AIC
#' @param mdl.class The name of the model's class
#' @param null.fixef Numeric vector containing the fixed effects of the null model.
#'        Only necessary for "poisson" family
#' @return A data frame with "Class", "Family", "Marginal", "Conditional", and "AIC" columns
.rsquared.glmm <- function(varF, varRand, varResid = NULL, varDisp = NULL, family, link,
                           mdl.aic, mdl.class, null.fixef = NULL){
  if(family == "gaussian"){
    # Only works with identity link
    if(link != "identity")
      family_link.stop(family, link)
    # Calculate marginal R-squared (fixed effects/total variance)
    Rm <- varF/(varF+varRand+varResid)
    # Calculate conditional R-squared (fixed effects+random effects/total variance)
    Rc <- (varF+varRand)/(varF+varRand+varResid)
  }
  else if(family == "binomial"){
    # Get the distribution-specific variance
    if(link == "logit")
      varDist <- (pi^2)/3
    else if(link == "probit")
      varDist <- 1
    else
      family_link.stop(family, link)
    # Calculate marginal R-squared
    Rm <- varF/(varF+varRand+varDist+varDisp)
    # Calculate conditional R-squared (fixed effects+random effects/total variance)
    Rc <- (varF+varRand)/(varF+varRand+varDist+varDisp)
  }
  else if(family == "poisson"){
    # Get the distribution-specific variance
    if(link == "log")
      varDist <- log(1+1/exp(null.fixef))
    else if(link == "sqrt")
      varDist <- 0.25
    else
      family_link.stop(family, link)
    # Calculate marginal R-squared
    Rm <- varF/(varF+varRand+varDist+varDisp)
    # Calculate conditional R-squared (fixed effects+random effects/total variance)
    Rc <- (varF+varRand)/(varF+varRand+varDist+varDisp)
  }
  else
    family_link.stop(family, link)
  # Bind R^2s into a matrix and return with AIC values
  data.frame(Class=mdl.class, Family = family, Link = link,
             Marginal=Rm, Conditional=Rc, AIC=mdl.aic)
}
 
#' stop execution if unable to calculate variance for a given family and link
family_link.stop <- function(family, link){
  stop(paste("Don't know how to calculate variance for",
             family, "family and", link, "link."))
}
