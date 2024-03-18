# package loading
library(haven)
library(tidyverse)
library(geepack)
library(lavaan)
library(WeightIt)
library(cobalt)

# Example for 1-period lagged effect --------------------------------------
#temporal relationship between systolic blood pressure and fasting blood glucose

# standardize the fast blood glucose:glu, systolic blood pressure:systolic
biomarker_dat <- mutate_at(biomarker_dat, 
                           .vars = c("glu", "systolic"),
                           .funs = function(x){(x-mean(x))/sd(x)})


# glu~systolic
dat_gs2 <- cbind(select(biomarker_dat, id, glu, systolic.2), fit_gs_2.trim[["weights"]], rep(2, nrow(biomarker_dat)))
colnames(dat_gs2) <- c("id", "glu", "systolic", "weight", "wave")
dat_gs1 <- cbind(select(biomarker_dat, id, glu.2, systolic.1), fit_gs_1.trim[["weights"]], rep(1, nrow(biomarker_dat)))
colnames(dat_gs1) <- c("id", "glu", "systolic", "weight", "wave")
dat_gee_gs <- rbind(dat_gs2,dat_gs1)
dat_gee_gs <- dat_gee_gs[order(dat_gee_gs$id, dat_gee_gs$wave), ]

fit_gs <- geeglm(glu ~ systolic, id = id, family = gaussian, data = dat_gee_gs,
                 weights = weight, corstr = "independence")
summary(fit_gs)

# systolic~glu
dat_sg2 <- cbind(select(biomarker_dat, id, systolic, glu.2), fit_sg_2.trim[["weights"]], rep(2, nrow(biomarker_dat)))
colnames(dat_sg2) <- c("id", "systolic", "glu", "weight", "wave")
dat_sg1 <- cbind(select(biomarker_dat, id, systolic.2, glu.1), fit_sg_1.trim[["weights"]], rep(1, nrow(biomarker_dat)))
colnames(dat_sg1) <- c("id", "systolic", "glu", "weight", "wave")
dat_gee_sg <- rbind(dat_sg2,dat_sg1)
dat_gee_sg <- dat_gee_sg[order(dat_gee_sg$id, dat_gee_sg$wave), ]

fit_sg <- geeglm(systolic ~ glu, id = id, family = gaussian, data = dat_gee_sg,
                 weights = weight, corstr = "independence")
summary(fit_sg)


# Example for 2-period lagged effect --------------------------------------
#temporal relationship between cognitive function and psychological well-being

# standardize the psychological well-being:personality, cognitive function: cognitive
person_cong_dat <- mutate_at(person_cong_dat, 
                             .vars = c("personality.1", "cognitive.1", 
                                       "personality.2", "cognitive.2",  
                                       "personality.3", "cognitive.3", 
                                       "personality.4", "cognitive.4"),
                             .funs = function(x){(x-mean(x))/sd(x)})

## weighting method selection
method <- c("glm", "gbm", "cbps", "npcbps","super")
fit <- vector("list", 5)
names(fit) <- method
for (i in 1:length(method)) {
  fit[[i]] <- vector("list", 2)
  names(fit[[i]]) c("pc", "cp")
}

formula_pc <- list(cognitive.2 ~ cognitive.1+personality.2+personality.1+
                     age.2+residenc.2+marriage.2+living.2+smoking.2+alcohol.2+diet.2+exercise.2+sex+ethnic+education,
                   cognitive.3 ~ cognitive.2+personality.3+personality.2+
                     age.3+residenc.3+marriage.3+living.3+smoking.3+alcohol.3+diet.3+exercise.3+sex+ethnic+education)
formula_cp <- list(personality.2 ~ personality.1+cognitive.2+cognitive.1+
                     age.2+residenc.2+marriage.2+living.2+smoking.2+alcohol.2+diet.2+exercise.2+sex+ethnic+education, 
                   personality.3 ~ personality.2+cognitive.3+cognitive.2+
                     age.3+residenc.3+marriage.3+living.3+smoking.3+alcohol.3+diet.3+exercise.3+sex+ethnic+education)

## 1.glm
fit[[1]][[1]] <- weightitMSM(formula_pc, data = person_cong_dat, method = "glm", stabilize = T)
fit[[1]][[2]] <- weightitMSM(formula_cp, data = person_cong_dat, method = "glm", stabilize = T)

## 2.gbm Propensity Score Weighting Using Generalized Boosted Models
fit[[2]][[1]] <- weightitMSM(formula_pc, data = person_cong_dat, method = "gbm", stabilize = T, criterion="p.mean")
fit[[2]][[2]] <- weightitMSM(formula_cp, data = person_cong_dat, method = "gbm", stabilize = T, criterion="p.mean")

## 3.cbps Covariate Balancing Propensity Score weighting
fit[[3]][[1]] <- weightitMSM(formula_pc, data = person_cong_dat, method = "cbps", over=F, stabilize = T)
fit[[3]][[2]] <- weightitMSM(formula_cp, data = person_cong_dat, method = "cbps", over=F, stabilize = T)

## 4.npcbps Non-parametric Covariate Balancing Propensity Score weighting
fit[[4]][[1]] <- weightitMSM(formula_pc, data = person_cong_dat, method = "npcbps", stabilize = T)
fit[[4]][[2]] <- weightitMSM(formula_cp, data = person_cong_dat, method = "npcbps", stabilize = T)

## 5.super Propensity Score Weighting Using SuperLearner
SL.library <- c("SL.glm", "SL.bartMachine", "SL.gam", "SL.glmnet", "SL.earth")
fit[[5]][[1]] <- weightitMSM(formula_pc,
                             data = person_cong_dat, method = "super", stabilize = T, cvControl = list(V = 20),
                             SL.library = SL.library)
fit[[5]][[2]] <- weightitMSM(formula_cp, 
                             data = person_cong_dat, method = "super", stabilize = T, cvControl = list(V = 20),
                             SL.library = SL.library)

#trim
fit.trim <- fit
for(i in 1:5){
  fit.trim[[i]][[1]] <- trim(fit[[i]][[1]], at=0.01, lower = T)
  fit.trim[[i]][[2]] <- trim(fit[[i]][[2]], at=0.01, lower = T)
}

balance <- summary <- vector("list", 5)
names(balance) <- names(summary) <- method
for (i in 1:length(method)) {
  balance[[i]] <- summary[[i]] <- vector("list", 2)
  names(balance[[i]]) <- names(summary[[i]]) <- c("pc", "cp")
}

for (i in 1:5) {
  for (j in 1:2) {
    balance[[i]][[j]] <- bal.tab(fit.trim[[i]][[j]])
  }
}

# estimating for the 5 weighting method
for (i in 1:length(method)) {
  # person~cong
  summary[[i]][[1]] <- summary(geeglm(personality.4 ~ cognitive.3+cognitive.2, 
                                      weights=fit.trim[[i]][[1]][["weights"]], 
                                      data =person_cong_dat, id= rownames(person_cong_dat)))
  # cong~person
  summary[[i]][[2]] <- summary(geeglm(cognitive.4 ~ personality.3+personality.2, 
                                      weights=fit.trim[[i]][[2]][["weights"]], 
                                      data =person_cong_dat, id= rownames(person_cong_dat)))
}
