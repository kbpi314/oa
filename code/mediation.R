table = read.csv('~/Desktop/clemente_lab/Projects/oa/inputs/df_med_R.tsv', sep = '\t')

independent = 'Lachno'
mediator = 'H_Nonhydroxylated_BA'
outcome = 'Pain'

independent = as.character(independent)
mediator = as.character(mediator)
outcome = as.character(outcome)
run_table <- table[!is.na(table[[independent]]) & table[[independent]] != "nan" &
                     !is.na(table[[mediator]]) & table[[mediator]] != "nan" &
                     !is.na(table[[outcome]]) & table[[outcome]] != "nan",
                   c(independent, mediator, outcome)]

                   
YX = summary(lm(run_table[[outcome]] ~ run_table[[independent]])) #, data = run_table))
MX = summary(lm(run_table[[mediator]] ~ run_table[[independent]]))#, data = run_table))
YM = summary(lm(run_table[[outcome]] ~ run_table[[mediator]]))#, data = run_table))
YMX = summary(lm(run_table[[outcome]] ~ run_table[[mediator]] + run_table[[independent]]))#, data = run_table))

# install.packages('mediation')
library(mediation)

model.m = lm(run_table[[mediator]] ~ run_table[[independent]])
model.y = lm(run_table[[outcome]] ~ run_table[[mediator]] + run_table[[independent]])

# M ~ X
form.m <- formula(paste("`", mediator, "` ~ `", independent, "`", sep=""))
model.m <- do.call(what="lm", list(formula=form.m, data=run_table))
# Y ~ M + X
form.y <- formula(paste("`", outcome, "` ~ `", mediator, "` + `", independent, "`", sep=""))
model.y <- do.call(what="lm", list(formula=form.y, data=run_table))


df <- data.frame()

for (i in seq(1:100)){
  med <- suppressMessages(mediation::mediate(model.m,
                                             model.y,
                                             boot=TRUE,
                                             treat=independent,
                                             mediator=mediator,
                                             sims=200))
  
  
  split = NA
  independent_sparsity <- nrow(run_table[run_table[[independent]] == 0,]) / nrow(run_table)
  mediator_sparsity <- nrow(run_table[run_table[[mediator]] == 0,]) / nrow(run_table)
  outcome_sparsity <- nrow(run_table[run_table[[outcome]] == 0,]) / nrow(run_table)
  
  results.df <- data.frame(X=independent,
                           M=mediator,
                           Y=outcome,
                           split_value=as.character(split),
                           total_effect=as.numeric(med$tau.coef),
                           total_effect_p=as.numeric(med$tau.p),
                           total_effect_ci_lower=as.numeric(med$tau.ci[1]),
                           total_effect_ci_upper=as.numeric(med$tau.ci[2]),
                           direct_effect=as.numeric(med$z0),
                           direct_effect_p=as.numeric(med$z0.p),
                           direct_effect_ci_lower=as.numeric(med$z0.ci[1]),
                           direct_effect_ci_upper=as.numeric(med$z0.ci[2]),
                           indirect_effect=as.numeric(med$d0),
                           indirect_effect_p=as.numeric(med$d0.p),
                           indirect_effect_ci_lower=as.numeric(med$d0.ci[1]),
                           indirect_effect_ci_upper=as.numeric(med$d0.ci[2]),
                           prop_mediated=as.numeric(med$n0),
                           prop_mediated_p=as.numeric(med$n0.p),
                           prop_mediated_ci_lower=as.numeric(med$n0.ci[1]),
                           prop_mediated_ci_upper=as.numeric(med$n0.ci[2]),
                           b2_effect=as.numeric(model.m$coefficients[2]),
                           b2_effect_p=as.numeric(summary(model.m)$coefficients[2,4]),
                           b2_effect_ci_lower=as.numeric(confint(model.m)[2,1]),
                           b2_effect_ci_upper=as.numeric(confint(model.m)[2,2]),
                           b3_effect=as.numeric(model.y$coefficients[3]),
                           b3_effect_p=as.numeric(summary(model.y)$coefficients[3,4]),
                           b3_effect_ci_lower=as.numeric(confint(model.y)[3,1]),
                           b3_effect_ci_upper=as.numeric(confint(model.y)[3,2]),
                           X_sparsity=as.numeric(independent_sparsity),
                           M_sparsity=as.numeric(mediator_sparsity),
                           Y_sparsity=as.numeric(outcome_sparsity),
                           n=as.numeric(nrow(run_table)))          

  df <- rbind(df, results.df)
}
write.csv(df,'~/Desktop/clemente_lab/Projects/oa/outputs/jobs20/df200.csv')

