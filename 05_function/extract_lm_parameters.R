
# extract q, r2 and coefficient from lm model


library(tidyverse)
stat_table <- function(glm.fit, qval_max, r2_min) {
  r2 <- glm.fit %>% map_df(~.$r.squared) %>% t()
  # retrieve p value for each fit (calculated based on F value)
  fval <- glm.fit %>% map_df(~.$fstatistic)
  pval <- map_df(fval, ~ pf(.[1], .[2], .[3], lower.tail = F)) %>% t()
  # calculate q value to correct p value
  qval <- p.adjust(pval, method = "fdr")
  # get a table containing gene ID as rowname, and r2, pval, qval
  stat_table_r2_q <-data.frame(r2, pval, qval)
  stat_table_r2_q <- add_column(stat_table_r2_q, ENSMUSG_ID= rownames(stat_table_r2_q), .before=1)
  #keep all the good fit genes (FDR < 0.05, r2 > 0.5) based on the predictor, then sort all the table based on q value (ascending)
  stat_table_r2_q <- stat_table_r2_q %>% filter(qval<qval_max, r2>r2_min) %>% arrange(desc(r2))
  #get all the good fit genes in stat_table_r2_q from the glm.fit, get all the "Estimate" and their "Pr(>|t|)" for downstream filtering
  fit.qr <- glm.fit[stat_table_r2_q$ENSMUSG_ID]
  # fit.qr[[1]] %>% coef() # take a look at the fit.qr
  fit.est <- fit.qr %>% map_df( ~ coef(.)[,"Estimate"]) %>%t()
  colnames(fit.est) <-fit.qr[[1]] %>% coef() %>% rownames() %>% paste("est", sep="_")
  fit.pval <- fit.qr %>% map_df(~ coef(.)[, "Pr(>|t|)"]) %>%t()
  colnames(fit.pval) <-fit.qr[[1]] %>% coef() %>% rownames() %>% paste("pval", sep="_")
  stat_table_all <- data.frame(stat_table_r2_q, fit.est, fit.pval)
  return(stat_table_all)
}
