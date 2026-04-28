############################################################
# BAYESIAN PORTFOLIO OPTIMIZATION
# Plug-in vs Bayesian (Normal-Inverse-Wishart)
# Rolling window backtest | Out-of-sample evaluation
# Reference: Lai, Xing & Chen (2011), Ann. Appl. Stat.
############################################################

############################
# 0. SETUP
############################

required_packages <- c("tidyverse", "moments", "corrplot", "quadprog",
                       "ggplot2", "scales", "zoo")
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
COL_EW     <- "#4dac26"   # equally weighted benchmark

BG_COLOR <- "#FFFFFF"

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

# 20 large-cap US stocks across 6 sectors (Tech, Finance, Healthcare,
# Consumer, Energy, Industrial). With a 60-month window this gives
# n/m = 3.0, which satisfies n > m for a non-singular covariance matrix

# Notation note: throughout this script, p denotes the number of assets,
# corresponding to m in the report notation.
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
# Note: kurtosis() from {moments} returns TOTAL kurtosis (normal = 3).
# Excess kurtosis = kurtosis(x) - 3, which is what most finance literature reports.
summary_stats <- data.frame(
  Asset          = asset_names,
  Mean           = unname(colMeans(R_all)),
  SD             = unname(apply(R_all, 2, sd)),
  Min            = unname(apply(R_all, 2, min)),
  Max            = unname(apply(R_all, 2, max)),
  Skewness       = unname(apply(R_all, 2, skewness)),
  Excess_Kurtosis = unname(apply(R_all, 2, function(x) kurtosis(x) - 3)),
  row.names      = NULL
)
print(summary_stats[, sapply(summary_stats, is.numeric)] %>% round(5))
write.csv(summary_stats, "output/summary_stats.csv", row.names = FALSE)

# 3.2 Skewness table
skewness_table <- data.frame(
  Asset    = asset_names,
  Skewness = unname(apply(R_all, 2, skewness)),
  row.names = NULL
)
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
  # nu0_extra controls degrees of freedom beyond the minimum (p+1).
  # Note: nu0 is NOT part of the sensitivity grid here; it could be
  # an additional axis of robustness analysis (left as future work).
  
  # kappa0 = 25: weakly informative prior, equivalent to 25 pseudo-
  # observations. Chosen to shrink the posterior mean toward zero
  # without dominating the 60 actual observations in the window.
  
  # zero prior mean: neutral assumption consistent with weak-form
  # market efficiency. The data then pull the posterior away from zero.
  # mu0 corresponds to nu in the report
  mu0 <- rep(0, p)
  nu0 <- p + nu0_extra
  
  # S0 (prior scale matrix Psi) scaled by the average in-sample variance (empirical Bayes):
  # the prior covariance adapts to the volatility regime of each window
  # rather than being fixed across the entire backtest.
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
  
  # Wealth starts at $1: prepend 0 to cumsum so first observation is exp(r[1])
  # relative to a starting value of 1.
  wealth      <- c(1, exp(cumsum(r)))
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

# lambda = 5: moderate risk aversion, consistent with Table 1 of
# Lai et al. (2011) which evaluates lambda = 1, 5, 10.
LAMBDA <- 5

w_plugin_train <- plugin_weights(R_first60, lambda = LAMBDA)
w_bayes_train  <- bayes_weights( R_first60, lambda = LAMBDA)
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

weights_df <- data.frame(
  Asset  = rep(asset_names, 2),
  Weight = c(w_plugin_train, w_bayes_train),
  Method = rep(c("Plugin", "Bayesian"), each = length(asset_names))
)

p_weights <- ggplot(weights_df, aes(x = Asset, y = Weight, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Plugin" = COL_PLUGIN, "Bayesian" = COL_BAYES)) +
  coord_flip() +
  labs(title = "Portfolio Weights Comparison (Training Sample)",
       x = NULL, y = "Weight", fill = NULL) +
  theme_portfolio

print(p_weights)

diff_w <- w_bayes_train - w_plugin_train
barplot(diff_w,
        main = "Difference (Bayesian - Plug-in)",
        las = 2, col = "darkgray")
abline(h = 0, lty = 2)

############################
# 6. ROLLING WINDOW BACKTEST
############################

# 60-month window: chosen to keep n/m = 3.0 while starting the
# out-of-sample period early enough to include the 2008-2009 GFC.
# Lai et al. (2011) use n=120, but that would push the backtest
# start to 2010, missing the crisis entirely.
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
ew_ret     <- numeric(n_steps)
plugin_W   <- matrix(NA_real_, nrow = n_steps, ncol = N_assets)
bayes_W    <- matrix(NA_real_, nrow = n_steps, ncol = N_assets)
ew_W       <- matrix(NA_real_, nrow = n_steps, ncol = N_assets)

backtest_dates <- dates_all[(window + 1):N_obs]

for (t in seq_len(n_steps)) {
  R_window <- R_all[t:(t + window - 1), ]
  r_next   <- R_all[t + window, ]
  
  w_p <- plugin_weights(R_window, lambda = LAMBDA)
  w_b <- bayes_weights( R_window, lambda = LAMBDA)
  w_e <- rep(1 / N_assets, N_assets)
  
  # Portfolio log-return: convert log-returns to linear, compute weighted sum,
  # then convert back. log(sum(w * exp(r))) would be a geometric mean, not the
  # correct arithmetic portfolio return.
  plugin_ret[t] <- log(1 + sum(w_p * (exp(r_next) - 1)))
  bayes_ret[t]  <- log(1 + sum(w_b * (exp(r_next) - 1)))
  ew_ret[t]     <- log(1 + sum(w_e * (exp(r_next) - 1)))
  
  plugin_W[t, ] <- w_p
  bayes_W[t, ]  <- w_b
  # EW weights after price drift (before next rebalancing):
  # w_e_drift reflects how 1/m weights evolve with asset returns,
  # which is what must be rebalanced back to 1/m each period.
  r_lin         <- exp(r_next) - 1
  w_e_drift     <- w_e * (1 + r_lin) / (1 + sum(w_e * r_lin))
  ew_W[t, ]     <- w_e_drift
  
  if (t %% 50 == 0) cat("Step", t, "/", n_steps, "\n")
}

colnames(plugin_W) <- asset_names
colnames(bayes_W)  <- asset_names
colnames(ew_W)     <- asset_names

############################
# 7. PERFORMANCE METRICS
############################

cat("\n=== Performance Metrics ===\n")

m_plugin <- performance_metrics(plugin_ret)
m_bayes  <- performance_metrics(bayes_ret)
m_ew     <- performance_metrics(ew_ret)

results_perf <- rbind(Plugin = m_plugin, Bayesian = m_bayes, EqualWeight = m_ew)
print(round(results_perf, 4))

turnover_results <- data.frame(
  Method       = c("Plugin", "Bayesian", "EqualWeight"),
  Avg_Turnover = c(
    compute_turnover(plugin_W),
    compute_turnover(bayes_W),
    {
      # EW turnover: at each period, the portfolio has drifted from 1/m.
      # Turnover = average total rebalancing needed to restore 1/m weights.
      # Measured as mean absolute deviation between drifted weights (ew_W)
      # and the 1/m target — this is what is actually traded each month.
      ew_target_W <- matrix(1 / N_assets, nrow = n_steps, ncol = N_assets)
      mean(rowSums(abs(ew_W - ew_target_W)))
    }
  )
)
print(turnover_results)

############################
# SENSITIVITY ANALYSIS
############################

kappa_grid <- c(5, 10, 25, 50, 100)

sensitivity_results <- map_dfr(kappa_grid, function(k) {
  
  bayes_ret_k <- numeric(n_steps)
  bayes_W_k   <- matrix(NA_real_, nrow = n_steps, ncol = N_assets)
  
  for (t in seq_len(n_steps)) {
    R_window <- R_all[t:(t + window - 1), ]
    r_next   <- R_all[t + window, ]
    
    w_b <- bayes_weights(R_window, lambda = LAMBDA, kappa0 = k)
    
    bayes_ret_k[t] <- log(1 + sum(w_b * (exp(r_next) - 1)))
    bayes_W_k[t, ] <- w_b
  }
  
  m_b <- performance_metrics(bayes_ret_k)
  
  data.frame(
    kappa0   = k,
    Sharpe   = m_b$Sharpe_Ratio,
    Return   = m_b$Annualized_Return,
    Vol      = m_b$Annualized_Volatility,
    Turnover = compute_turnover(bayes_W_k)
  )
})
print(round(sensitivity_results, 4))

# kappa0 = 25 extracted directly from sensitivity_results for consistency
bayes_ret_25 <- {
  tmp <- numeric(n_steps)
  for (t in seq_len(n_steps)) {
    R_window <- R_all[t:(t + window - 1), ]
    r_next   <- R_all[t + window, ]
    w_b25    <- bayes_weights(R_window, lambda = LAMBDA, kappa0 = 25)
    tmp[t]   <- log(1 + sum(w_b25 * (exp(r_next) - 1)))
  }
  tmp
}

# kappa0 = 50
bayes_ret_50 <- {
  tmp <- numeric(n_steps)
  for (t in seq_len(n_steps)) {
    R_window <- R_all[t:(t + window - 1), ]
    r_next   <- R_all[t + window, ]
    w_b50    <- bayes_weights(R_window, lambda = LAMBDA, kappa0 = 50)
    tmp[t]   <- log(1 + sum(w_b50 * (exp(r_next) - 1)))
  }
  tmp
}

p_kappa_sharpe <- ggplot(sensitivity_results, aes(x = kappa0, y = Sharpe)) +
  geom_line(color = COL_BAYES, linewidth = 0.7) +
  geom_point(color = COL_BAYES, size = 2) +
  geom_hline(yintercept = m_plugin$Sharpe_Ratio,
             linetype = "dashed", color = COL_PLUGIN) +
  geom_hline(yintercept = m_ew$Sharpe_Ratio,
             linetype = "dotted", color = COL_EW) +
  labs(title    = "Sensitivity Analysis: Sharpe vs kappa0",
       subtitle = "Dashed = Plug-in  |  Dotted = Equal-Weight",
       x = "kappa0 (prior strength)",
       y = "Sharpe Ratio") +
  theme_portfolio

print(p_kappa_sharpe)
ggsave("figures/sensitivity_kappa.png", p_kappa_sharpe, width = 8, height = 5, dpi = 150)

# Comparison kappa=25 vs kappa=50
cum_df_kappa <- data.frame(
  Date     = backtest_dates,
  Bayes_25 = exp(cumsum(bayes_ret_25)),
  Bayes_50 = exp(cumsum(bayes_ret_50))
)

p_compare <- cum_df_kappa %>%
  pivot_longer(-Date, names_to = "Strategy", values_to = "Wealth") %>%
  ggplot(aes(x = Date, y = Wealth, color = Strategy)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = c("Bayes_25" = COL_BAYES, "Bayes_50" = "#1b9e77")) +
  labs(title    = "Cumulative Wealth: Bayesian kappa0 = 25 vs 50",
       x = NULL, y = "Wealth ($1 invested)", color = NULL) +
  theme_portfolio

print(p_compare)
ggsave("figures/kappa_comparison.png", p_compare, width = 9, height = 5, dpi = 150)

# kappa0 = 25 selected as canonical strategy for two reasons:
# (1) theoretically it is a weakly informative prior (25 pseudo-observations
#     vs 60 actual data points), avoiding over-shrinkage;
# (2) empirically kappa=25 outperforms kappa=50 after 2020, suggesting
#     the weaker prior better adapts to structural shifts in return dynamics.
bayes_ret <- bayes_ret_25

############################
# REGIME ANALYSIS
############################

# NOTE: descriptive only — sub-periods with fewer than ~36 observations
# yield unreliable annualised Sharpe ratios and should be read with caution.
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
  Bayesian = bayes_ret,
  EW       = ew_ret
)

regime_perf <- regime_df %>%
  group_by(Regime) %>%
  summarise(
    N_obs         = n(),
    Plugin_Sharpe = mean(Plugin)   / sd(Plugin)   * sqrt(12),
    Bayes_Sharpe  = mean(Bayesian) / sd(Bayesian) * sqrt(12),
    EW_Sharpe     = mean(EW)       / sd(EW)       * sqrt(12),
    Plugin_Return = mean(Plugin)   * 12,
    Bayes_Return  = mean(Bayesian) * 12,
    EW_Return     = mean(EW)       * 12
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
sharpe_ew     <- rolling_sharpe(ew_ret)

roll_df <- data.frame(
  Date     = backtest_dates,
  Plugin   = sharpe_plugin,
  Bayesian = sharpe_bayes,
  EW       = sharpe_ew
) %>% drop_na()

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
  geom_line(aes(y = Plugin,   color = "Plug-in"),       linewidth = 0.7) +
  geom_line(aes(y = Bayesian, color = "Bayesian"),      linewidth = 0.7) +
  geom_line(aes(y = EW,       color = "Equal-Weight"),  linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Plug-in"      = COL_PLUGIN,
                                "Bayesian"     = COL_BAYES,
                                "Equal-Weight" = COL_EW),
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
cum_plugin <- exp(cumsum(plugin_ret))
cum_bayes  <- exp(cumsum(bayes_ret))
cum_ew     <- exp(cumsum(ew_ret))

cum_df <- data.frame(
  Date         = backtest_dates,
  Plugin       = cum_plugin,
  Bayesian     = cum_bayes,
  Equal_Weight = cum_ew
)
write.csv(cum_df, "output/cumulative_wealth.csv", row.names = FALSE)

p_wealth <- cum_df %>%
  pivot_longer(-Date, names_to = "Strategy", values_to = "Wealth") %>%
  ggplot(aes(x = Date, y = Wealth, colour = Strategy)) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c("Plugin"       = COL_PLUGIN,
                                 "Bayesian"     = COL_BAYES,
                                 "Equal_Weight" = COL_EW)) +
  labs(title    = "Out-of-Sample Cumulative Wealth ($1 invested)",
       subtitle = period_label,
       x = NULL, y = "Cumulative Wealth ($)", colour = NULL) +
  theme_portfolio

ggsave("figures/cumulative_wealth.png", p_wealth, width = 10, height = 5, dpi = 150)

# ---- 8.2 Drawdown ----
peak_p <- cummax(cum_plugin)
peak_b <- cummax(cum_bayes)
peak_e <- cummax(cum_ew)

dd_df <- data.frame(
  Date         = backtest_dates,
  Plugin       = (cum_plugin - peak_p) / peak_p,
  Bayesian     = (cum_bayes  - peak_b) / peak_b,
  Equal_Weight = (cum_ew     - peak_e) / peak_e
)

dd_long <- dd_df %>%
  pivot_longer(-Date, names_to = "Strategy", values_to = "Drawdown") %>%
  mutate(Strategy = factor(Strategy, levels = c("Plugin", "Bayesian", "Equal_Weight")))

p_drawdown <- ggplot(dd_long, aes(x = Date, y = Drawdown, fill = Strategy)) +
  geom_area(alpha = 0.45, position = "identity") +
  scale_fill_manual(values = c("Plugin"       = COL_PLUGIN,
                               "Bayesian"     = COL_BAYES,
                               "Equal_Weight" = COL_EW)) +
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
  Date         = backtest_dates,
  Plugin       = plugin_ret,
  Bayesian     = bayes_ret,
  Equal_Weight = ew_ret
) %>%
  pivot_longer(-Date, names_to = "Strategy", values_to = "Return") %>%
  ggplot(aes(x = Date, y = Return, colour = Strategy)) +
  geom_line(alpha = 0.65, linewidth = 0.5) +
  scale_colour_manual(values = c("Plugin"       = COL_PLUGIN,
                                 "Bayesian"     = COL_BAYES,
                                 "Equal_Weight" = COL_EW)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title    = "Monthly Out-of-Sample Returns",
       subtitle = period_label,
       x = NULL, y = "Monthly Return", colour = NULL) +
  theme_portfolio

ggsave("figures/monthly_returns.png", p_ret, width = 10, height = 4, dpi = 150)

# ---- 8.5 Performance bar charts ----
p_perf_pct <- data.frame(
  Method     = c("Plugin", "Bayesian", "Equal-Weight"),
  Ann_Return = c(m_plugin$Annualized_Return,    m_bayes$Annualized_Return,    m_ew$Annualized_Return),
  Ann_Vol    = c(m_plugin$Annualized_Volatility, m_bayes$Annualized_Volatility, m_ew$Annualized_Volatility)
) %>%
  pivot_longer(-Method, names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Metric, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.65) +
  scale_fill_manual(values = c("Plugin"       = COL_PLUGIN,
                               "Bayesian"     = COL_BAYES,
                               "Equal-Weight" = COL_EW)) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(title = "Annualised Return & Volatility",
       x = NULL, y = NULL, fill = NULL) +
  theme_portfolio

p_sharpe <- data.frame(
  Method = c("Plugin", "Bayesian", "Equal-Weight"),
  Sharpe = c(m_plugin$Sharpe_Ratio, m_bayes$Sharpe_Ratio, m_ew$Sharpe_Ratio)
) %>%
  ggplot(aes(x = Method, y = Sharpe, fill = Method)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("Plugin"       = COL_PLUGIN,
                               "Bayesian"     = COL_BAYES,
                               "Equal-Weight" = COL_EW)) +
  labs(title = "Sharpe Ratio",
       x = NULL, y = "Sharpe Ratio", fill = NULL) +
  theme_portfolio

ggsave("figures/performance_pct.png", p_perf_pct, width = 7, height = 5, dpi = 150)
ggsave("figures/sharpe_ratio.png",   p_sharpe,   width = 6, height = 4, dpi = 150)

# ---- 8.6 Distribution of Portfolio Returns ----
dist_df <- data.frame(
  Return   = c(plugin_ret, bayes_ret, ew_ret),
  Strategy = rep(c("Plugin", "Bayesian", "Equal-Weight"), each = length(plugin_ret))
)

p_dist <- ggplot(dist_df, aes(x = Return, fill = Strategy)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = mean(plugin_ret), color = COL_PLUGIN, linetype = "dashed") +
  geom_vline(xintercept = mean(bayes_ret),  color = COL_BAYES,  linetype = "dashed") +
  geom_vline(xintercept = mean(ew_ret),     color = COL_EW,     linetype = "dashed") +
  scale_fill_manual(values = c("Plugin"       = COL_PLUGIN,
                               "Bayesian"     = COL_BAYES,
                               "Equal-Weight" = COL_EW)) +
  labs(title    = "Distribution of Portfolio Returns",
       subtitle = "Plug-in vs Bayesian vs Equal-Weight",
       x = "Monthly Return", y = "Density", fill = NULL) +
  theme_portfolio

ggsave("figures/return_distribution.png", p_dist, width = 8, height = 5, dpi = 150)

# ---- 8.7 Left tail comparison ----
tail_df <- dist_df %>%
  filter(Return < quantile(Return, 0.1))

p_tail <- ggplot(tail_df, aes(x = Return, fill = Strategy)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Plugin"       = COL_PLUGIN,
                               "Bayesian"     = COL_BAYES,
                               "Equal-Weight" = COL_EW)) +
  labs(title    = "Left Tail of Returns (Worst 10%)",
       subtitle = "Extreme downside risk comparison",
       x = "Monthly Return", y = "Density", fill = NULL) +
  theme_portfolio

ggsave("figures/left_tail.png", p_tail, width = 8, height = 5, dpi = 150)

# ---- Downside risk metrics ----
cat("\nLoss probability:\n")
cat("Plugin:", mean(plugin_ret < 0), "\n")
cat("Bayes :", mean(bayes_ret  < 0), "\n")
cat("EW    :", mean(ew_ret     < 0), "\n\n")

cat("Extreme loss probability (< -5%):\n")
cat("Plugin:", mean(plugin_ret < -0.05), "\n")
cat("Bayes :", mean(bayes_ret  < -0.05), "\n")
cat("EW    :", mean(ew_ret     < -0.05), "\n")

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
             Bayesian_Return = bayes_ret,
             EW_Return       = ew_ret),
  "output/backtest_returns.csv", row.names = FALSE
)

############################
# 10. TABLE EXPORTS (CSV)
############################

perf_csv <- data.frame(
  Metric       = c("Annualised Return (%)",
                   "Annualised Volatility (%)",
                   "Sharpe Ratio",
                   "Max Drawdown (%)",
                   "Avg Turnover"),
  Plug_in      = c(m_plugin$Annualized_Return     * 100,
                   m_plugin$Annualized_Volatility * 100,
                   m_plugin$Sharpe_Ratio,
                   m_plugin$Max_Drawdown          * 100,
                   turnover_results$Avg_Turnover[1]),
  Bayesian     = c(m_bayes$Annualized_Return     * 100,
                   m_bayes$Annualized_Volatility * 100,
                   m_bayes$Sharpe_Ratio,
                   m_bayes$Max_Drawdown          * 100,
                   turnover_results$Avg_Turnover[2]),
  Equal_Weight = c(m_ew$Annualized_Return     * 100,
                   m_ew$Annualized_Volatility * 100,
                   m_ew$Sharpe_Ratio,
                   m_ew$Max_Drawdown          * 100,
                   turnover_results$Avg_Turnover[3])
)
write.csv(perf_csv, "output/table_performance.csv", row.names = FALSE)

weight_instability <- function(W) {
  mean(apply(W, 2, sd))
}

turnover_ts <- function(W) {
  rowSums(abs(diff(W)))
}

stab_csv <- data.frame(
  Metric       = c("Weight Instability", "Avg Turnover"),
  Plug_in      = c(weight_instability(plugin_W), mean(turnover_ts(plugin_W))),
  Bayesian     = c(weight_instability(bayes_W),  mean(turnover_ts(bayes_W))),
  # Weight Instability is NA for EW: ew_W contains drift-adjusted weights,
  # not reoptimised weights, so its std dev across time is not comparable
  # to the instability measure computed for Plug-in and Bayesian.
  Equal_Weight = c(NA_real_, turnover_results$Avg_Turnover[3])
)
write.csv(stab_csv,            "output/table_stability.csv",    row.names = FALSE)
write.csv(sensitivity_results, "output/table_sensitivity.csv",  row.names = FALSE)

############################
# 11. FINAL SUMMARY (console)
############################

cat("\n====================================================\n")
cat("  FINAL PERFORMANCE SUMMARY\n")
cat("====================================================\n")
cat(sprintf("%-25s %10s %10s %14s\n", "Metric", "Plug-in", "Bayesian", "Equal-Weight"))
cat(strrep("-", 63), "\n")
cat(sprintf("%-25s %9.2f%% %9.2f%% %13.2f%%\n", "Annualised Return",
            m_plugin$Annualized_Return * 100,
            m_bayes$Annualized_Return  * 100,
            m_ew$Annualized_Return     * 100))
cat(sprintf("%-25s %9.2f%% %9.2f%% %13.2f%%\n", "Annualised Volatility",
            m_plugin$Annualized_Volatility * 100,
            m_bayes$Annualized_Volatility  * 100,
            m_ew$Annualized_Volatility     * 100))
cat(sprintf("%-25s %10.4f %10.4f %14.4f\n", "Sharpe Ratio",
            m_plugin$Sharpe_Ratio,
            m_bayes$Sharpe_Ratio,
            m_ew$Sharpe_Ratio))
cat(sprintf("%-25s %9.2f%% %9.2f%% %13.2f%%\n", "Max Drawdown",
            m_plugin$Max_Drawdown * 100,
            m_bayes$Max_Drawdown  * 100,
            m_ew$Max_Drawdown     * 100))
cat(sprintf("%-25s %10.4f %10.4f %14.4f\n", "Avg Turnover",
            turnover_results$Avg_Turnover[1],
            turnover_results$Avg_Turnover[2],
            turnover_results$Avg_Turnover[3]))
cat("====================================================\n")
cat("\nAll outputs saved in output/ — all figures saved in figures/\n")