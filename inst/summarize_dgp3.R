devtools::load_all()
library(ggplot2)
library(dplyr)

## summarize results of run_dgp3.R

# load saved results
experiment_result <- readRDS("script/result/dgp3_2sample_sumstat.rds")
list2env(experiment_result, .GlobalEnv)
beta_star <- c(beta1, beta2, 0, 0, 0)
alpha <- 0.05

# extract results
df_lst_est <- lapply(df_lst, \(ll) ll[[1]])
df_lst_res <- lapply(df_lst, \(ll) ll[[2]])

# combine results
df_res <- do.call(rbind, df_lst_res)
df_res$n <- factor(as.numeric(df_res$n))
df_res$method <- factor(df_res$method, levels = c("Q-subset-sumstat", 
                                                  "IV-Lasso-sumstat-refit", 
                                                  "IV-sumstat"))
df_res$valid_ass <- df_res$valid_ass1 & df_res$valid_ass2
table(df_res$valid_ass1)/(length(unique(df_res$method))*length(n_grid))
table(df_res$valid_ass2)/(length(unique(df_res$method))*length(n_grid))

df_res_long <- df_res[df_res$valid_ass,] %>%
  tidyr::pivot_longer(cols = 1:2, names_to = "name", values_to = "value")

p1 <- ggplot(df_res_long %>%
               dplyr::filter(method %in% c("Q-subset-sumstat",
                                           "IV-Lasso-sumstat-refit",
                                           "IV-sumstat")),
             aes(x = n, y = value, 
                 group = interaction(name, method, n),
                 color = method)) +
  ggbeeswarm::geom_quasirandom(dodge.width = 0.75, size = 0.5, alpha = 0.75) +
  geom_boxplot(width = 0.5, position = position_dodge(0.75), outlier.alpha = 0,
               alpha = 0) + 
  facet_grid(name ~.) +
  scale_y_log10() + 
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        plot.margin = margin(0, 0, 0, 0)) +
  scale_color_manual(values = c("#5E4FA2", "#74C476", "#FDAE61"),
                     labels = c("Q-subset-sumstat" = "spaceTSIV-L0",
                                "IV-Lasso-sumstat-refit" = "spaceTSIV-L1",
                                "IV-sumstat" = "TSIV"),
                     name = "estimator") +
  theme_bw() +
  ggplot2::labs(y = element_blank()) + 
  theme(legend.position = "top",
        strip.background = element_blank(),
        legend.box.margin = margin(0, 0, -1, 0), 
        legend.spacing.y = unit(-1, "pt"),
        legend.margin = margin(-4, 0, -4, 0),
        text = element_text(size = 13)) +
  ggplot2::guides(color = guide_legend(nrow = 1, byrow = FALSE))

df_est <- do.call(rbind, df_lst_est) %>%
  dplyr::filter(df_res$valid_ass) %>%
  dplyr::filter(method %in% c("Q-subset-sumstat",
                              "IV-Lasso-sumstat-refit",
                              "IV-sumstat")) %>%
  tidyr::pivot_longer(1:5, names_to = "beta_hat", values_to = "est") %>%
  dplyr::mutate(n = factor(n),
                method = factor(method, levels = c("Q-subset-sumstat",
                                                   "IV-Lasso-sumstat-refit",
                                                   "IV-sumstat"))) %>%
  dplyr::group_by(method, n) %>%
  dplyr::mutate(iter = rep(reps, each = length(beta_star)),
                star = rep(beta_star, length(reps)))

p3 <- ggplot(df_est, 
             aes(x = n, y = est, color = beta_hat)) + 
  geom_boxplot(width = 0.5, position = position_dodge(0.75), outlier.alpha = 0) + 
  ggbeeswarm::geom_quasirandom(dodge.width = 0.75, size = 0.5, alpha = 0.75) +
  geom_hline(linetype = 2, aes(yintercept = star, color = beta_hat)) + 
  facet_grid(method~., labeller = labeller(method = c("Q-subset-sumstat" = "spaceTSIV-L0",
                                                      "IV-Lasso-sumstat-refit" = "spaceTSIV-L1",
                                                      "IV-sumstat" = "TSIV"))) +
  xlab("sample size") + 
  ylab("Estimates") +
  scale_color_manual(values = colorspace::qualitative_hcl(length(unique(df_est$beta_hat)),
                                                          alpha = 0.75)) +
  coord_cartesian(ylim = c(-1, 3)) +
  labs(subtitle = "Note: some outliers are excluded but box stats are calculated using all results.") + 
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 9), 
        strip.background = element_blank(),
        strip.text = element_text(size = 9))

ggsave(paste0("script/result/dgp3_2sample.pdf"), 
       plot = p1, width = 5, height = 3.5)

ggsave(paste0("script/result/dgp3_2sample_estimates.pdf"), 
       plot = p3, width = 6, height = 7)

# count (percentage) of each sparsity
count_sparsity <- function(estimated_sparsity, divide_by = 100) {
  res <- table(estimated_sparsity)
  tibble::tibble(sparsity = names(res),
                 count = as.numeric(res)/divide_by)
}
df_out <- df_res[df_res$valid_ass,] %>%
  dplyr::group_by(method, n) %>%
  dplyr::summarise(count_sparsity(estimated_sparsity, divide_by = length(reps))) %>%
  tidyr::pivot_wider(id_cols = 1:2, 
                     names_from = sparsity, 
                     names_prefix = "sparsity = ",
                     names_sort = TRUE, 
                     values_from = count, 
                     values_fill = 0)

data.table::fwrite(x = df_out, 
                   paste0("script/result/dgp3_2sample.csv"))

# confusion matrix
confusion_matrix <- df_est %>%
  dplyr::filter(method != "IV-Lasso-sumstat-refit") %>%
  dplyr::group_by(n, method, beta_hat) %>%
  dplyr::summarise(tpr = sum((est != 0) & (star != 0)) / length(reps),
                   fpr = sum((est != 0) & (star == 0)) / length(reps),
                   fnr = sum((est == 0) & (star != 0)) / length(reps),
                   tnr = sum((est == 0) & (star == 0)) / length(reps))

# calculate jaccard 
res_jaccard <- df_est %>%
  dplyr::filter(!method %in% c("IV-sumstat")) %>%
  dplyr::group_by(method, n, iter) %>%
  dplyr::summarise(jaccard = jaccard_dist(which(est != 0), which(star != 0))) %>%
  dplyr::group_by(method, n) %>%
  dplyr::summarise(est = mean(jaccard),
                   lwr = min(mean(jaccard) - sd(jaccard), 1),
                   upr = min(mean(jaccard) + sd(jaccard), 1))
res_jaccard$name <- "jaccard"
res_jaccard$alpha <- NA

# pct of correct sparsity: 1 if the sparsity is correct and 0 otherwise 
res_sparsity <- df_est %>% 
  dplyr::group_by(method, n, iter) %>%
  dplyr::summarise(correct_sparsity = as.integer(sum(est != 0) == sum(star != 0))) %>%
  dplyr::mutate(correct_sparsity = if_else(is.na(correct_sparsity), 
                                           0L, correct_sparsity)) %>%
  dplyr::group_by(method, n) %>%
  dplyr::summarise(binom_test(sum(correct_sparsity), n()))
res_sparsity$name <- "correct_sparsity_pct"
res_sparsity$alpha <- 0.05 * 100
res_sparsity$est <- res_sparsity$est * 100
res_sparsity$lwr <- res_sparsity$lwr * 100
res_sparsity$upr <- res_sparsity$upr * 100

# calculate fpr
res_fpr <- df_est %>%
  dplyr::filter(!method %in% c("IV-sumstat"),
                star == 0) %>%
  dplyr::group_by(method, n, iter, beta_hat) %>%
  dplyr::mutate(fp = as.numeric((est != 0) & (star == 0))) %>%
  dplyr::group_by(method, n) %>%
  dplyr::summarise(binom_test(sum(fp), n(), alpha = 0.05)) 
res_fpr$name <- "fpr"
res_fpr$alpha <- 0.05

# calculate tpr
res_tpr <- df_est %>%
  dplyr::filter(!method %in% c("IV-sumstat"),
                star != 0) %>%
  dplyr::group_by(method, n, iter, beta_hat) %>%
  dplyr::mutate(tp = as.numeric((est != 0) & (star != 0))) %>%
  dplyr::group_by(method, n) %>%
  dplyr::summarise(binom_test(sum(tp), n(), alpha = 0.05)) 
res_tpr$name <- "tpr"
res_tpr$alpha <- 0.05 * 100
res_tpr$est <- res_tpr$est * 100
res_tpr$lwr <- res_tpr$lwr * 100
res_tpr$upr <- res_tpr$upr * 100

select_res <- rbind(res_jaccard, res_sparsity, res_tpr) %>%
  dplyr::filter(method %in% c("Q-subset-sumstat", "IV-Lasso-sumstat-refit")) %>%
  dplyr::mutate(name = factor(name, 
                              levels = c("jaccard", "correct_sparsity_pct", "tpr")))

p4 <- select_res %>%
  ggplot2::ggplot(aes(x = n, 
                      y = est, 
                      group = interaction(name, method),
                      color = method, 
                      fill = method)) +
  ggplot2::geom_line(position = position_dodge(0.2)) +
  ggplot2::geom_hline(aes(yintercept = alpha), 
                      color = 'red', 
                      size = 0.5, 
                      alpha = 0.8,
                      linetype = 2) +
  ggplot2::geom_errorbar(aes(ymin = lwr, ymax = upr), width = .15,
                         position = position_dodge(0.2)) +
  ggplot2::facet_grid(name ~., scale = "free_y", 
                      labeller = as_labeller(c("jaccard" = "Jaccard similarity", 
                                               "correct_sparsity_pct" = "correct size (%)",
                                               "tpr" = "tpr (%)"))) + 
  ggplot2::labs(y = element_blank()) + 
  ggplot2::scale_color_manual(values = c("#5E4FA2", "#74C476", "#FDAE61"),
                              labels = c("Q-subset-sumstat" = "spaceTSIV-L0",
                                         "IV-Lasso-sumstat-refit" = "spaceTSIV-L1"),
                              name = "estimator") + 
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "top",
                 strip.background = element_blank(),
                 legend.box.margin = margin(0, 0, -1, 0), 
                 legend.spacing.y = unit(-1, "pt"),
                 legend.margin = margin(-4, 0, -4, 0),
                 text = element_text(size = 13),
                 strip.text = element_text(size = 7))

p <- ggpubr::ggarrange(p1, p4, nrow = 1, common.legend = FALSE)

ggsave(paste0("script/result/dgp3_2sample_combined.pdf"), 
       plot = p, width = 10, height = 3.5)
