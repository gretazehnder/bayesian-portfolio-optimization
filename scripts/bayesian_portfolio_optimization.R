############################################################
# BAYESIAN PORTFOLIO OPTIMIZATION
# Plug-in vs Bayesian (Normal-Inverse-Wishart)
# Rolling window backtest | Out-of-sample evaluation
# Reference: Lai, Xing & Chen (2011), Ann. Appl. Stat.
############################################################

############################
# 0. SETUP
############################

required_packages <- c("tidyverse", "moments", "corrplot", "quadprog", "ggplot2", "scales", "zoo")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(tidyverse)
library(moments)
library(corrplot)
library(quadprog)
library(ggplot2)
library(scales)
library(zoo)

dir.create("output",  showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

COL_PLUGIN <- "#d73027"
COL_BAYES  <- "#2166ac"

# ---- Chart background theme ----
BG_COLOR   <- "#F5F5F0"

theme_portfolio <- theme_minimal(base_size = 13) +
  theme(
    legend.position   = "top",
    plot.title        = element_text(face = "bold"),
    plot.subtitle     = element_text(colour = "grey40"),
    plot.background   = element_rect(fill = BG_COLOR, color = NA),
    panel.background  = element_rect(fill = BG_COLOR, color = NA),
    legend.background = element_rect(fill = BG_COLOR, color = NA)
  )

############################
# 1. LOAD DATA
############################

stocks      <- read.csv("data/stocks.csv", stringsAsFactors = FALSE)
stocks$Date <- as.Date(stocks$Date)
stocks      <- stocks %>% arrange(Date)

cat("\n=== Data ===\n")
cat("Rows:   ", nrow(stocks), "\n")
cat("Assets: ", ncol(stocks) - 1, "\n")
cat("Period: ", as.character(min(stocks$Date)),
    "to", as.character(max(stocks$Date)), "\n")

# ---- NA check ----
na_check <- colSums(is.na(stocks %>% select(-Date)))
na_check <- na_check[na_check > 0]
if (length(na_check) == 0) {
  cat("NA check: no missing values found.\n\n")
} else {
  cat("NAs found — applying forward fill:\n")
  print(na_check)
  cat("\n")
  stocks <- stocks %>%
    mutate(across(where(is.numeric), ~ zoo::na.locf(., na.rm = FALSE)))
}

############################
# 2. LOG RETURNS
############################

returns <- stocks %>%
  mutate(across(where(is.numeric), ~ log(. / lag(.)))) %>%
  drop_na()

R_all       <- returns %>% select(-Date) %>% as.matrix()
dates_all   <- returns$Date
N_assets    <- ncol(R_all)
asset_names <- colnames(R_all)
N_obs       <- nrow(R_all)

cat("Monthly obs: ", N_obs, "\n")
cat("Assets:      ", N_assets, "\n\n")

############################
# 3. EDA  (full sample)
############################

cat("=== EDA ===\n")

# 3.1 Summary statistics
summary_stats <- data.frame(
  Asset    = asset_names,
  Mean     = unname(colMeans(R_all)),
  SD       = unname(apply(R_all, 2, sd)),
  Min      = unname(apply(R_all, 2, min)),
  Max      = unname(apply(R_all, 2, max)),
  Skewness = unname(apply(R_all, 2, skewness)),
  Kurtosis = unname(apply(R_all, 2, kurtosis)),
  row.names = NULL
)
print(summary_stats[, sapply(summary_stats, is.numeric)] %>% round(5))
write.csv(summary_stats, "output/summary_stats.csv", row.names = FALSE)

# 3.2 Skewness table
skewness_table <- data.frame(
  Asset     = asset_names,
  Skewness  = unname(apply(R_all, 2, skewness)),
  row.names = NULL
)
skewness_table
write.csv(skewness_table, "output/skewness_table.csv", row.names = FALSE)

# 3.3 Shapiro-Wilk normality test
# H0: normal distribution  |  reject H0 if p-value < 0.05
shapiro_table <- data.frame(
  Asset           = asset_names,
  Shapiro_p_value = unname(apply(R_all, 2, function(x) shapiro.test(x)$p.value)),
  row.names       = NULL
)
cat("Normal assets (Shapiro p > 0.05):",
    sum(shapiro_table$Shapiro_p_value > 0.05), "/", N_assets, "\n")
shapiro_table
write.csv(shapiro_table, "output/shapiro_table.csv", row.names = FALSE)

# 3.4 Correlation & covariance (full sample)
correlation_matrix <- cor(R_all)
covariance_matrix  <- cov(R_all)
write.csv(as.data.frame(correlation_matrix), "output/correlation_matrix.csv", row.names = TRUE)
write.csv(as.data.frame(covariance_matrix),  "output/covariance_matrix.csv",  row.names = TRUE)

png("figures/correlation_matrix.png", width = 1000, height = 800)
par(bg = BG_COLOR)
corrplot(correlation_matrix, method = "color", type = "upper", tl.cex = 0.7,
         title = "Correlation Matrix — Full Sample (2000-2025)",
         mar   = c(0, 0, 2, 0),
         col   = colorRampPalette(c(COL_PLUGIN, "white", COL_BAYES))(200))
dev.off()

# 3.5 Correlation matrix — first training window (60 months)
R_first60 <- R_all[1:60, ]
png("figures/correlation_matrix_training.png", width = 1000, height = 800)
par(bg = BG_COLOR)
corrplot(cor(R_first60), method = "color", type = "upper", tl.cex = 0.7,
         title = "Correlation Matrix — First Training Window (60 months)",
         mar   = c(0, 0, 2, 0),
         col   = colorRampPalette(c(COL_PLUGIN, "white", COL_BAYES))(200))
dev.off()

############################
# 4. FUNCTIONS
############################

# ------------------------------------------------------------------
# 4.1  Plug-in: classical Markowitz, long-only, fully invested
# ------------------------------------------------------------------
plugin_weights <- function(R_window, lambda = 5) {
  p         <- ncol(R_window)
  mu_hat    <- colMeans(R_window)
  Sigma_hat <- cov(R_window) + diag(1e-6, p)
  
  Dmat <- 2 * lambda * Sigma_hat
  dvec <- mu_hat
  Amat <- cbind(rep(1, p), diag(p))
  bvec <- c(1, rep(0, p))
  meq  <- 1
  
  sol <- tryCatch(
    solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = meq),
    error = function(e) NULL
  )
  if (is.null(sol)) return(rep(1 / p, p))
  w <- pmax(sol$solution, 0)
  w / sum(w)
}

# ------------------------------------------------------------------
# 4.2  Bayesian NIW: Normal-Inverse-Wishart conjugate prior
# ------------------------------------------------------------------
bayes_weights <- function(R_window, lambda = 5,
                          kappa0    = 25,
                          nu0_extra = 20) {
  X  <- as.matrix(R_window)
  n  <- nrow(X)
  p  <- ncol(X)
  
  xbar     <- colMeans(X)
  S_sample <- cov(X)
  
  # ---- Prior hyperparameters ----
  mu0 <- rep(0, p)
  nu0 <- p + nu0_extra
  
  prior_var <- mean(diag(S_sample))
  S0 <- prior_var * diag(p) * (nu0 - p - 1)
  
  # ---- Posterior update ----
  kappa_n <- kappa0 + n
  nu_n    <- nu0 + n
  mu_n    <- (kappa0 * mu0 + n * xbar) / kappa_n
  diff_mu <- matrix(xbar - mu0, ncol = 1)
  S_n     <- S0 +
    (n - 1) * S_sample +
    (kappa0 * n / kappa_n) * (diff_mu %*% t(diff_mu))
  
  # ---- Posterior predictive moments ----
  scale_factor <- (kappa_n + 1) / (kappa_n * (nu_n - p - 1))
  Sigma_pred   <- scale_factor * S_n + diag(1e-6, p)
  
  Dmat <- 2 * lambda * Sigma_pred
  dvec <- mu_n
  Amat <- cbind(rep(1, p), diag(p))
  bvec <- c(1, rep(0, p))
  meq  <- 1
  
  sol <- tryCatch(
    solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = meq),
    error = function(e) NULL
  )
  
  if (is.null(sol)) return(rep(1 / p, p))
  w <- pmax(sol$solution, 0)
  w / sum(w)
}

# ------------------------------------------------------------------
# 4.3  Performance metrics (MONTHLY data)
# ------------------------------------------------------------------
performance_metrics <- function(r) {
  ann_return <- mean(r) * 12
  ann_vol    <- sd(r)   * sqrt(12)
  sharpe     <- ifelse(ann_vol > 0, ann_return / ann_vol, NA)
  wealth      <- cumprod(1 + r)
  running_max <- cummax(wealth)
  drawdown    <- (wealth - running_max) / running_max
  max_dd      <- min(drawdown)
  
  data.frame(
    Annualized_Return     = ann_return,
    Annualized_Volatility = ann_vol,
    Sharpe_Ratio          = sharpe,
    Max_Drawdown          = max_dd
  )
}

# ------------------------------------------------------------------
# 4.4  Turnover
# ------------------------------------------------------------------
compute_turnover <- function(W) {
  if (nrow(W) < 2) return(NA)
  mean(rowSums(abs(diff(W))))
}

############################
# 5. TRAINING SAMPLE WEIGHTS
############################

w_plugin_train <- plugin_weights(R_first60, lambda = 5)
w_bayes_train  <- bayes_weights( R_first60, lambda = 5)
names(w_plugin_train) <- asset_names
names(w_bayes_train)  <- asset_names

cat("\n=== Training Weights (first 60 months) ===\n")
cat("Plug-in:  "); print(round(w_plugin_train, 4))
cat("Bayesian: "); print(round(w_bayes_train,  4))

png("figures/plugin_weights_training.png", width = 1200, height = 700, bg = BG_COLOR)
par(bg = BG_COLOR)
barplot(w_plugin_train, main = "Plug-in Weights (Training Sample)",
        las = 2, col = COL_PLUGIN,
        ylim = c(0, max(w_plugin_train) * 1.3), ylab = "Weight")
dev.off()

png("figures/bayesian_weights_training.png", width = 1200, height = 700, bg = BG_COLOR)
par(bg = BG_COLOR)
barplot(w_bayes_train, main = "Bayesian Weights (Training Sample)",
        las = 2, col = COL_BAYES,
        ylim = c(0, max(w_bayes_train) * 1.3), ylab = "Weight")
dev.off()

graphics.off()
par(mfrow = c(1, 2))

barplot(w_plugin_train,
        main = "Plug-in Weights",
        las = 2, col = COL_PLUGIN,
        ylim = c(0, max(c(w_plugin_train, w_bayes_train)) * 1.3))

barplot(w_bayes_train,
        main = "Bayesian Weights",
        las = 2, col = COL_BAYES,
        ylim = c(0, max(c(w_plugin_train, w_bayes_train)) * 1.3))

par(mfrow = c(1, 1))

diff_w <- w_bayes_train - w_plugin_train
barplot(diff_w,
        main = "Difference (Bayesian - Plug-in)",
        las = 2, col = "darkgray")
abline(h = 0, lty = 2)

############################
# 6. ROLLING WINDOW BACKTEST
############################

window  <- 60
n_steps <- N_obs - window

cat("\n=== Rolling Window Backtest ===\n")
cat("Window:           ", window, "months\n")
cat("Rebalancing steps:", n_steps, "\n")
cat("Out-of-sample:    ",
    as.character(dates_all[window + 1]), "to",
    as.character(dates_all[N_obs]), "\n\n")

plugin_ret <- numeric(n_steps)
bayes_ret  <- numeric(n_steps)
plugin_W   <- matrix(NA_real_, nrow = n_steps, ncol = N_assets)
bayes_W    <- matrix(NA_real_, nrow = n_steps, ncol = N_assets)

backtest_dates <- dates_all[(window + 1):N_obs]

for (t in seq_len(n_steps)) {
  R_window <- R_all[t:(t + window - 1), ]
  r_next   <- R_all[t + window, ]
  
  w_p <- plugin_weights(R_window, lambda = 5)
  w_b <- bayes_weights( R_window, lambda = 5)
  
  plugin_ret[t] <- sum(w_p * r_next)
  bayes_ret[t]  <- sum(w_b * r_next)
  plugin_W[t, ] <- w_p
  bayes_W[t, ]  <- w_b
  
  if (t %% 50 == 0) cat("Step", t, "/", n_steps, "\n")
}

colnames(plugin_W) <- asset_names
colnames(bayes_W)  <- asset_names

############################
# 7. PERFORMANCE METRICS
############################

cat("\n=== Performance Metrics ===\n")

m_plugin <- performance_metrics(plugin_ret)
m_bayes  <- performance_metrics(bayes_ret)

results_perf <- rbind(Plugin = m_plugin, Bayesian = m_bayes)
print(round(results_perf, 4))

turnover_results <- data.frame(
  Method       = c("Plugin", "Bayesian"),
  Avg_Turnover = c(compute_turnover(plugin_W),
                   compute_turnover(bayes_W))
)
print(turnover_results)

############################
# REGIME ANALYSIS
############################

regime <- case_when(
  backtest_dates <= as.Date("2007-12-31") ~ "Expansion",
  backtest_dates <= as.Date("2009-12-31") ~ "Global Financial Crisis",
  backtest_dates <  as.Date("2020-02-01") ~ "Low Vol / Expansion",
  backtest_dates <  as.Date("2021-12-31") ~ "COVID Crisis & Recovery",
  backtest_dates <  as.Date("2023-01-01") ~ "Inflation / Rate Shock",
  TRUE                                    ~ "Recent Stabilization"
)

regime_df <- data.frame(
  Date     = backtest_dates,
  Regime   = regime,
  Plugin   = plugin_ret,
  Bayesian = bayes_ret
)

regime_perf <- regime_df %>%
  group_by(Regime) %>%
  summarise(
    Plugin_Sharpe = mean(Plugin)   / sd(Plugin)   * sqrt(12),
    Bayes_Sharpe  = mean(Bayesian) / sd(Bayesian) * sqrt(12),
    Plugin_Return = mean(Plugin)   * 12,
    Bayes_Return  = mean(Bayesian) * 12
  )

print(regime_perf)

############################
# ROLLING SHARPE RATIO
############################

rolling_sharpe <- function(r, window = 24) {
  rollapply(r, width = window, FUN = function(x) {
    if (sd(x) == 0) return(NA)
    mean(x) / sd(x) * sqrt(12)
  }, fill = NA, align = "right")
}

sharpe_plugin <- rolling_sharpe(plugin_ret)
sharpe_bayes  <- rolling_sharpe(bayes_ret)

roll_df <- data.frame(
  Date     = backtest_dates,
  Plugin   = sharpe_plugin,
  Bayesian = sharpe_bayes
) %>% drop_na()

# Crisis regime bands
regime_bands <- data.frame(
  xmin  = as.Date(c("2007-12-01", "2020-02-01", "2022-01-01")),
  xmax  = as.Date(c("2009-12-31", "2021-12-31", "2022-12-31")),
  label = c("GFC", "COVID", "Inflation"),
  fill  = c("#f4a9a8", "#ffe0b2", "#e8d5f5")
)

p_rolling_sharpe <- ggplot(roll_df, aes(x = Date)) +
  geom_rect(data = regime_bands,
            aes(xmin = xmin, xmax = xmax,
                ymin = -Inf, ymax = Inf, fill = label),
            alpha = 0.4, inherit.aes = FALSE) +
  scale_fill_manual(values = setNames(regime_bands$fill, regime_bands$label),
                    name = "Regime") +
  geom_line(aes(y = Plugin,   color = "Plug-in"),  linewidth = 0.9) +
  geom_line(aes(y = Bayesian, color = "Bayesian"), linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Plug-in"  = COL_PLUGIN,
                                "Bayesian" = COL_BAYES),
                     name = "Strategy") +
  labs(title    = "Rolling Sharpe Ratio (24-month window)",
       subtitle = "Out-of-sample performance over time",
       y = "Sharpe Ratio", x = NULL) +
  theme_portfolio

ggsave("figures/rolling_sharpe.png", p_rolling_sharpe, width = 11, height = 5, dpi = 150)
print(p_rolling_sharpe)

############################
# 8. FIGURES
############################

cat("\n=== Generating Figures ===\n")

period_label <- paste0(
  format(backtest_dates[1],       "%b %Y"), " - ",
  format(tail(backtest_dates, 1), "%b %Y")
)

# ---- 8.1 Cumulative Wealth ----
cum_plugin <- cumprod(1 + plugin_ret)
cum_bayes  <- cumprod(1 + bayes_ret)

cum_df <- data.frame(
  Date     = backtest_dates,
  Plugin   = cum_plugin,
  Bayesian = cum_bayes
)
write.csv(cum_df, "output/cumulative_wealth.csv", row.names = FALSE)

p_wealth <- cum_df %>%
  pivot_longer(-Date, names_to = "Strategy", values_to = "Wealth") %>%
  ggplot(aes(x = Date, y = Wealth, colour = Strategy)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("Plugin"   = COL_PLUGIN,
                                 "Bayesian" = COL_BAYES)) +
  labs(title    = "Out-of-Sample Cumulative Wealth ($1 invested)",
       subtitle = period_label,
       x = NULL, y = "Cumulative Wealth ($)", colour = NULL) +
  theme_portfolio

ggsave("figures/cumulative_wealth.png", p_wealth, width = 10, height = 5, dpi = 150)

# ---- 8.2 Drawdown (Bayesian plotted first so Plug-in is on top) ----
peak_p <- cummax(cum_plugin)
peak_b <- cummax(cum_bayes)

dd_df <- data.frame(
  Date     = backtest_dates,
  Plugin   = (cum_plugin - peak_p) / peak_p,
  Bayesian = (cum_bayes  - peak_b) / peak_b
)

# Pivot with Bayesian first so it renders behind Plug-in
dd_long <- dd_df %>%
  pivot_longer(-Date, names_to = "Strategy", values_to = "Drawdown") %>%
  mutate(Strategy = factor(Strategy, levels = c("Plugin", "Bayesian")))

p_drawdown <- ggplot(dd_long, aes(x = Date, y = Drawdown, fill = Strategy)) +
  geom_area(alpha = 0.50, position = "identity") +
  scale_fill_manual(values = c("Plugin"   = COL_PLUGIN,
                               "Bayesian" = COL_BAYES)) +
  scale_y_continuous(
    labels = function(x) paste0(round(x * 100, 0), "%"),
    breaks = seq(-0.4, 0, by = 0.1)
  ) +
  labs(title    = "Drawdown",
       subtitle = period_label,
       x = NULL, y = "Drawdown", fill = NULL) +
  theme_portfolio

ggsave("figures/drawdown.png", p_drawdown, width = 10, height = 4, dpi = 150)

# ---- 8.3 Weight instability ----
p_wvol <- data.frame(
  Asset    = asset_names,
  Plugin   = unname(apply(plugin_W, 2, sd)),
  Bayesian = unname(apply(bayes_W,  2, sd)),
  row.names = NULL
) %>%
  pivot_longer(-Asset, names_to = "Method", values_to = "SD_Weight") %>%
  ggplot(aes(x = reorder(Asset, SD_Weight), y = SD_Weight, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Plugin" = COL_PLUGIN, "Bayesian" = COL_BAYES)) +
  coord_flip() +
  labs(title    = "Weight Instability: Std Dev of Allocation Over Time",
       subtitle = "Higher = more unstable rebalancing",
       x = NULL, y = "Std Dev of Weight", fill = NULL) +
  theme_portfolio

ggsave("figures/weight_instability.png", p_wvol, width = 9, height = 7, dpi = 150)

# ---- 8.4 Monthly returns ----
p_ret <- data.frame(
  Date     = backtest_dates,
  Plugin   = plugin_ret,
  Bayesian = bayes_ret
) %>%
  pivot_longer(-Date, names_to = "Strategy", values_to = "Return") %>%
  ggplot(aes(x = Date, y = Return, colour = Strategy)) +
  geom_line(alpha = 0.65, linewidth = 0.5) +
  scale_colour_manual(values = c("Plugin"   = COL_PLUGIN,
                                 "Bayesian" = COL_BAYES)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title    = "Monthly Out-of-Sample Returns",
       subtitle = period_label,
       x = NULL, y = "Monthly Return", colour = NULL) +
  theme_portfolio

ggsave("figures/monthly_returns.png", p_ret, width = 10, height = 4, dpi = 150)

# ---- 8.5 Performance bar chart ----
p_perf <- data.frame(
  Method     = c("Plugin", "Bayesian"),
  Ann_Return = c(m_plugin$Annualized_Return,    m_bayes$Annualized_Return),
  Ann_Vol    = c(m_plugin$Annualized_Volatility, m_bayes$Annualized_Volatility),
  Sharpe     = c(m_plugin$Sharpe_Ratio,          m_bayes$Sharpe_Ratio)
) %>%
  pivot_longer(-Method, names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Metric, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.65) +
  scale_fill_manual(values = c("Plugin"   = COL_PLUGIN,
                               "Bayesian" = COL_BAYES)) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(title = "Annualised Performance Metrics",
       x = NULL, y = NULL, fill = NULL) +
  theme_portfolio

ggsave("figures/performance_summary.png", p_perf, width = 7, height = 5, dpi = 150)

# ---- 8.6 Distribution of Portfolio Returns ----
dist_df <- data.frame(
  Return   = c(plugin_ret, bayes_ret),
  Strategy = rep(c("Plugin", "Bayesian"), each = length(plugin_ret))
)

p_dist <- ggplot(dist_df, aes(x = Return, fill = Strategy)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = mean(plugin_ret), color = COL_PLUGIN, linetype = "dashed") +
  geom_vline(xintercept = mean(bayes_ret),  color = COL_BAYES,  linetype = "dashed") +
  scale_fill_manual(values = c("Plugin"   = COL_PLUGIN,
                               "Bayesian" = COL_BAYES)) +
  labs(title    = "Distribution of Portfolio Returns",
       subtitle = "Plug-in vs Bayesian",
       x = "Monthly Return", y = "Density", fill = NULL) +
  theme_portfolio

ggsave("figures/return_distribution.png", p_dist, width = 8, height = 5, dpi = 150)

# ---- 8.7 Left tail comparison ----
tail_df <- dist_df %>%
  filter(Return < quantile(Return, 0.1))

p_tail <- ggplot(tail_df, aes(x = Return, fill = Strategy)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Plugin"   = COL_PLUGIN,
                               "Bayesian" = COL_BAYES)) +
  labs(title    = "Left Tail of Returns (Worst 10%)",
       subtitle = "Extreme downside risk comparison",
       x = "Monthly Return", y = "Density", fill = NULL) +
  theme_portfolio

ggsave("figures/left_tail.png", p_tail, width = 8, height = 5, dpi = 150)

# ---- Downside risk metrics ----
cat("\nLoss probability:\n")
cat("Plugin:", mean(plugin_ret < 0), "\n")
cat("Bayes :", mean(bayes_ret  < 0), "\n\n")

cat("Extreme loss probability (< -5%):\n")
cat("Plugin:", mean(plugin_ret < -0.05), "\n")
cat("Bayes :", mean(bayes_ret  < -0.05), "\n")

############################
# 9. SAVE OUTPUTS
############################

write.csv(summary_stats,    "output/summary_stats.csv",       row.names = FALSE)
write.csv(skewness_table,   "output/skewness_table.csv",      row.names = FALSE)
write.csv(shapiro_table,    "output/shapiro_table.csv",       row.names = FALSE)
write.csv(results_perf,     "output/performance_results.csv", row.names = TRUE)
write.csv(turnover_results, "output/turnover_results.csv",    row.names = FALSE)
write.csv(
  data.frame(Date            = backtest_dates,
             Plugin_Return   = plugin_ret,
             Bayesian_Return = bayes_ret),
  "output/backtest_returns.csv", row.names = FALSE
)

############################
# 10. FINAL SUMMARY
############################

cat("\n====================================================\n")
cat("  FINAL PERFORMANCE SUMMARY\n")
cat("====================================================\n")
cat(sprintf("%-25s %10s %10s\n", "Metric", "Plug-in", "Bayesian"))
cat(strrep("-", 47), "\n")
cat(sprintf("%-25s %9.2f%% %9.2f%%\n", "Annualised Return",
            m_plugin$Annualized_Return * 100,
            m_bayes$Annualized_Return  * 100))
cat(sprintf("%-25s %9.2f%% %9.2f%%\n", "Annualised Volatility",
            m_plugin$Annualized_Volatility * 100,
            m_bayes$Annualized_Volatility  * 100))
cat(sprintf("%-25s %10.4f %10.4f\n", "Sharpe Ratio",
            m_plugin$Sharpe_Ratio,
            m_bayes$Sharpe_Ratio))
cat(sprintf("%-25s %9.2f%% %9.2f%%\n", "Max Drawdown",
            m_plugin$Max_Drawdown * 100,
            m_bayes$Max_Drawdown  * 100))
cat(sprintf("%-25s %10.4f %10.4f\n", "Avg Turnover",
            turnover_results$Avg_Turnover[1],
            turnover_results$Avg_Turnover[2]))
cat("====================================================\n")
cat("\nOutputs -> output/\nFigures  -> figures/\nDone.\n")

############################
# 11. STABILITY METRICS
############################

weight_instability <- function(W) {
  mean(apply(W, 1, function(w) sum((w - mean(w))^2)))
}

turnover_ts <- function(W) {
  rowSums(abs(diff(W)))
}

cat("\n=== STABILITY METRICS ===\n")
cat("Weight Instability (Plugin): ", weight_instability(plugin_W), "\n")
cat("Weight Instability (Bayes):  ", weight_instability(bayes_W),  "\n")
cat("Avg Turnover (Plugin):       ", mean(turnover_ts(plugin_W)),  "\n")
cat("Avg Turnover (Bayes):        ", mean(turnover_ts(bayes_W)),   "\n")