# Title; Propensity Score Matching Analysis Template
# Author: Mark Freeman II
# Date: December 2020

# The goal of this r file is to provide a template for propensity score matching (PSM). PSM
# is an observational studies method that is useful for accounting for observed covariates in
# a dataset that you did not have control of and or for secondary analyses.

# highly enncourage reading the article "Matching methods for causal inference: A review 
# and a look forward" by Elizabeth Stuart for more information regarding propensity score
# matching: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2943670/

# process to conduct a propensity score matching analysis:
# step 0: import packages, data, and clean data
# step 1: create distance measures (e.g. propensity score)
# step 2: match observations using distance measure
# step 3: assess quality of match through covariate balance
# step 4: outcomes analysis

##############################
# step 0: importing packages #
##############################

library(tidyverse) # data cleaning
library(forcats) # working with categorical vars
library(tableone) # easy table 1 (i.e. summary stats for medical research)
library(kableExtra) # beautiful R markdown tables (HTML format)
library(optmatch) # package for propensity score matching
library(cobalt) # extension to psm packages for checking covariate balance with beautiful charts
library(stargazer) # cool package to share results of different regression models in one table

##########################
# step 0: importing data #
##########################

# file path and name
file_path <- '~/your/file/path/'
file_name <- 'file_name'

# strongly suggest using col_types = cols_only() argument to bring in specific values as desired data type
your_df <- read_csv(
  paste(file_path, file_name, sep=''),
  col_types = cols_only(
    'outcome_var' = col_integer(),
    'treatment_var' = col_character(),
    'covariates_1' = col_integer(),
    'covariates_2' = col_character(),
    'covariates_n' = col_guess()
    )
  )

#########################
# step 0: cleaning data #
#########################

# this process is dependent on dataset; follow CRISP-DM method

#########################################
# step 1: calculating propensity scores #
#########################################

# note that there are more than just propensity score for distance measures (e.g. mahalanobis distance)

# "treatment": treatment_var
# covariates: covariates_1, covariates_2, covariates_n
propensity_glm <- glm(treatment_var ~ covariates_1 + covariates_2 + covariates_n, data=your_df, family=binomial)

# propesnity score is successful if there is no violation of strongly ignorable treatment assignment (SITA)
# If violation of SITA, then remove a covariate(s) until error is resolved 9should be informed by domain knowledge
# and experimental design)

####################
# step 2: matching #
####################

# matching feels more like an art than a science in this step
# you will have to try multiple methods until you get a successful match
# highly encourage starting simple then moving to more involved methodes, as this impacts approach for
# outcomes analysis

# below are some different examples

# pair match is one of the simplest matching methods, and is preferred if match quality and experimental
# design allows
# pair match, 1:k ratio, k=3
match_result_k3 <- pairmatch(
  match_on(propensity_glm),
  data=your_df,
  controls=3
  )

# full match is very robust and often leads to great matches, but then requires extra methods for outcome analysis
# full match, many:many ratio, min.controls=2
match_result_full <- fullmatch(
  match_on(propensity_glm),
  data=your_df,
  min.controls=2
  )

# exact match within a full match is very useful for when a specifc covariate is very important for
# experimental design
# full match, many:many ratio, exact match on specific covariate
match_result_full_exact <- fullmatch(
  match_on(propensity_glm),
  data=your_df,
  within=exactMatch(
    treatment_var ~ covariates_2,
    data=your_df
    )
  )

###################################
# step 3: assess quality of match #
###################################

# quality of match is assessed via the standardized mean difference (SMD), where a SMD of <0.1 for all
# covariates is considered an excellent match

# using bal.tab from cobalt package
cov_balance <- bal.tab(
  match_result_k3,
  treatment_var ~ covariates_1 + covariates_2 + covariates_n,
  data=your_df,
  un=TRUE
  )

print(cov_balance)

# plotting covariate balance (i.e. absolute standardized mean difference) via cobalt package
love.plot(cov_balance,
          stat = 'mean.diffs',
          threshold = .1, 
          var.order = 'unadjusted',
          abs = TRUE,
          line = TRUE,
          limits = c(0, 1))

#############################
# step 4: outcomes analysis #
#############################

# very important note: if you use 1:k matching, then linear regression will suffice for
# outcomes analysis. if you used full matching, then you need to account for the matched group 
# differences through methods such as inverse probability weighting

# extract matched pairs from matching model
# subset data to only matched observations
matched_df <- your_df
matched_df$matched_pairs <- match_result_k3
matched_df_final <- matched_df %>%
  drop_na(matched_pairs)

# models comparing impact of matching

# simple OLS
unadj_ols <- lm(outcome_var ~ treatment_var, data=your_df)
summary(unadj_ols)

# multiple linear regression
unadj_lm <- lm(outcome_var ~ treatment_var + covariates_1 + covariates_2 + covariates_n, data=your_df)
summary(unadj_lm)

# matched data
matched_lm <- lm(outcome_var ~ treatment_var + covariates_1 + covariates_2 + covariates_n, data=matched_df_final)
summary(matched_lm)

# compare all models with stargazer
stargazer(
  unadj_ols,unadj_lm,matched_lm,
  intercept.bottom = FALSE,
  covariate.labels = c(
    'Intercept',
    'Treatment',
    'Covariate 1',
    'Covariate 2',
    'Covariate n'
    ),
  dep.var.caption  = 'Dependent Variable: Outcome Variable',
  dep.var.labels   = '',
  column.labels = c('OLS - No Match', 'Multiple Reg - No Match', 'Multiple Reg - Match'),
  type = 'text',
  model.numbers = FALSE,
  single.row = TRUE
  )

# matched_lm diagnostic plots
layout(matrix(c(1,2,3,4),2,2)) 
plot(matched_lm)

###################
# end of analysis #
###################
