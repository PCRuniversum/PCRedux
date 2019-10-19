library(qpcR)
library(PCRedux)

library(dplyr)
library(forcats)
library(reshape2)

library(mlr)

library(ggplot2)

if("devtools" %in% rownames(installed.packages()) == FALSE) {
  library(patchwork)
}

if("gbm" %in% rownames(installed.packages()) == FALSE) {
  install.packages("gbm")
}

if("patchwork" %in% rownames(installed.packages())) {
  library(patchwork)
} else {
  devtools::install_github("thomasp85/patchwork")
  "patchwork" %in% rownames(installed.packages())
  library(patchwork)
}
#####

curves <- 245

# visualizng qPCR curves
p1 <- data.frame(htPCR[, 1L:curves]) %>% 
  melt(id.vars = "Cycles") %>%
  ggplot(aes(x = Cycles, y = value, color = variable)) +
  geom_line() +
  theme_bw() +
  xlab("Cycle") +
  ylab("Raw fluorescence") +
  ggtitle("A)") + #qPCR curves
  theme(legend.position = "none")

# obtaining decisions
dec_htPCR <- read.csv(system.file("decision_res_htPCR.csv", package = "PCRedux"))
dec <- unlist(lapply(1L:(curves -1), function(i) {
  decision_modus(dec_htPCR[i, 2:8])
}))

# calculating encu parameters
# Please be patient, this step will take some time
res <- encu(htPCR[, 1L:curves])

# merging into one dataset
dat <- cbind(res, decision = factor(c("ambiguous", "negative", "positive")[dec], 
                                    levels = c("positive", "ambiguous", "negative"))) %>%
  select(eff, loglin_slope, minRFU, init2, decision) %>%
  filter(!is.na(dec))

# visualizing
p2 <- ggplot(data = dat %>%
               mutate(id = rownames(dat)) %>%
               melt(id.vars = c("id", "decision")), 
             aes(x = decision, y = value)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~ variable, scales = "free_y") +
  scale_y_continuous("Value") +
  scale_x_discrete("Assessment") +
  ggtitle("B)") # Separation of types of curves by encu() parameters

# modelling
tsk <- makeClassifTask("pcr_classif", data = dat, target = "decision")

mdls <- list()
mdls[[1]] <- makeLearner("classif.ranger", predict.type = "prob")
mdls[[2]] <- makeLearner("classif.ksvm", predict.type = "prob")
mdls[[3]] <- makeLearner("classif.lda", predict.type = "prob")
mdls[[4]] <- makeLearner("classif.gbm", predict.type = "prob")
mdls[[5]] <- makeLearner("classif.multinom", predict.type = "prob")
mdls[[6]] <- makeLearner("classif.glmnet", predict.type = "prob")

set.seed(4732)
results <- do.call(rbind, lapply(mdls, function(mdl) {
  res <- resample(mdl, tsk, cv10, measures = list(mmce, multiclass.au1u))
  cbind(model = res[["learner.id"]], res[["measures.test"]])
}))
results[["model"]] <- fct_recode(results[["model"]],
                                 `ranger::ranger` = "classif.ranger", 
                                 `kernlab::ksvm` = "classif.ksvm",
                                 `MASS::lda` = "classif.lda",
                                 `gbm::gbm` = "classif.gbm",
                                 `nnet::multinom` = "classif.multinom",
                                 `glmnet::glmnet` = "classif.glmnet")

p3 <- ggplot(data = results, aes(x = model, y = multiclass.au1u)) +
  geom_point() + 
  geom_errorbar(data = results %>% 
                  group_by(model) %>% 
                  summarise(auc = median(multiclass.au1u)), 
                aes(x = model, ymin = auc, ymax = auc), 
                inherit.aes = FALSE, color = "#FC5E61",
                width = 0.5) +
  theme_bw() + 
  xlab("Model") +
  ylab("AUC one vs all mean result") +
  ggtitle("C)") # Results of crossvalidating models trained on encu() parameters

cairo_ps("figure1.eps", width = 10, height = 5)
(p1 + p2) / p3
dev.off()
