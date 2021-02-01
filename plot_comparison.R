rm(list = ls()); gc(reset = TRUE)
# plot_path <- "yourpath"; setwd(plot_path)
library(reshape2)
library(ggplot2)

u_method_type <- c("REG", "REG-M",
                   "DRL", "DRL-M", "MIS", "MIS-M", "COPE")
method_type <- factor(rep(u_method_type, 2),
                      levels = u_method_type)

assessment_type <- rep(c("Bias", "MSE"), each = length(u_method_type))

result_summary1 <- function(reg_1, reg_2, dr, drl, mis, trl, trim = 0.00) {
  root_experiment_rep <- sqrt(ncol(reg_1))
  
  bias_tab <- rbind(apply(reg_1, 1, mean, trim = trim),
                    apply(reg_2, 1, mean, trim = trim),
                    apply(dr, 1, mean, trim = trim),
                    apply(drl, 1, mean, trim = trim),
                    apply(mis, 1, mean, trim = trim),
                    apply(trl, 1, mean, trim = trim))
  
  sd_tab <- rbind(apply(reg_1, 1, sd),
                  apply(reg_2, 1, sd),
                  apply(dr, 1, sd),
                  apply(drl, 1, sd),
                  apply(mis, 1, sd),
                  apply(trl, 1, sd))
  sd_tab <- sd_tab / root_experiment_rep
  
  mse <- function(x, trim = 0) {
    if (trim > 0) {
      q1 <- quantile(x, trim)
      q2 <- quantile(x, 1 - trim)
      x <- x[x >= q1 & x <= q2]
    }
    mean(x^2)
  }
  
  mse_tab <- rbind(apply(reg_1, 1, mse, trim = trim),
                   apply(reg_2, 1, mse, trim = trim),
                   apply(dr, 1, mse, trim = trim),
                   apply(drl, 1, mse, trim = trim),
                   apply(mis, 1, mse, trim = trim),
                   apply(trl, 1, mse, trim = trim))
  
  mssd_tab <- rbind(apply(reg_1^2, 1, sd),
                    apply(reg_2^2, 1, sd),
                    apply(dr^2, 1, sd),
                    apply(drl^2, 1, sd),
                    apply(mis^2, 1, sd),
                    apply(trl^2, 1, sd))
  mssd_tab <- mssd_tab / root_experiment_rep
  
  list(bias_tab, sd_tab, mse_tab, mssd_tab)
}

result_summary2 <- function(reg_1, reg_1_wm, dr, drl, mis, mis_wm, trl, trim = 0.00) {
  root_experiment_rep <- sqrt(ncol(reg_1))
  
  bias_tab <- rbind(apply(abs(reg_1), 1, mean, trim = trim),
                    apply(abs(reg_1_wm), 1, mean, trim = trim),
                    apply(abs(dr), 1, mean, trim = trim),
                    apply(abs(drl), 1, mean, trim = trim),
                    apply(abs(mis), 1, mean, trim = trim),
                    apply(abs(mis_wm), 1, mean, trim = trim),
                    apply(abs(trl), 1, mean, trim = trim))
  
  sd_tab <- rbind(apply(abs(reg_1), 1, sd),
                  apply(abs(reg_1_wm), 1, sd),
                  apply(abs(dr), 1, sd),
                  apply(abs(drl), 1, sd),
                  apply(abs(mis), 1, sd),
                  apply(abs(mis_wm), 1, sd),
                  apply(abs(trl), 1, sd))
  sd_tab <- sd_tab / root_experiment_rep
  
  mse <- function(x, trim = 0) {
    if (trim > 0) {
      q1 <- quantile(x, trim)
      q2 <- quantile(x, 1 - trim)
      x <- x[x >= q1 & x <= q2]
    }
    mean(x^2)
  }
  
  mse_tab <- rbind(apply(reg_1, 1, mse, trim = trim),
                   apply(reg_1_wm, 1, mse, trim = trim),
                   apply(dr, 1, mse, trim = trim),
                   apply(drl, 1, mse, trim = trim),
                   apply(mis, 1, mse, trim = trim),
                   apply(mis_wm, 1, mse, trim = trim),
                   apply(trl, 1, mse, trim = trim))
  
  mssd_tab <- rbind(apply(reg_1^2, 1, sd),
                    apply(reg_1_wm^2, 1, sd),
                    apply(dr^2, 1, sd),
                    apply(drl^2, 1, sd),
                    apply(mis^2, 1, sd),
                    apply(mis_wm^2, 1, sd),
                    apply(trl^2, 1, sd))
  mssd_tab <- mssd_tab / root_experiment_rep
  
  list(bias_tab, sd_tab, mse_tab, mssd_tab)
}

trajectory_num_vec <- c(60, 100, 140, 180, 220)
trajectory_num_vec_len <- length(trajectory_num_vec)

reg <- read.csv("trajectory-reg.csv", header = FALSE)
regwm <- read.csv("trajectory-regwm.csv", header = FALSE)
drl <- read.csv("trajectory-drl.csv", header = FALSE)
drlwm <- read.csv("trajectory-drlwm.csv", header = FALSE)
mis <- read.csv("trajectory-is.csv", header = FALSE)
miswm <- read.csv("trajectory-iswm.csv", header = FALSE)
trl <- read.csv("trajectory-trl.csv", header = FALSE)
res_trajectory <- result_summary2(reg, regwm, drl, drlwm, mis, miswm, trl)

reg <- read.csv("time-reg.csv", header = FALSE)
regwm <- read.csv("time-regwm.csv", header = FALSE)
drl <- read.csv("time-drl.csv", header = FALSE)
drlwm <- read.csv("time-drlwm.csv", header = FALSE)
mis <- read.csv("time-is.csv", header = FALSE)
miswm <- read.csv("time-iswm.csv", header = FALSE)
trl <- read.csv("time-trl.csv", header = FALSE)
res_time <- result_summary2(reg, regwm, drl, drlwm, mis, miswm, trl)

process_result <- function(mean_tab, sd_tab, row_name, col_name) {
  factored_row_name <- factor(row_name, levels = row_name)
  
  pdat1 <- as.data.frame(mean_tab)
  colnames(pdat1) <- col_name
  pdat1[["Method"]] <- factored_row_name
  pdat1 <- melt(pdat1, c("Method"))
  
  pdat2 <- as.data.frame(sd_tab)
  colnames(pdat2) <- col_name
  pdat2[["Method"]] <- factored_row_name
  pdat2 <- melt(pdat2, c("Method"))
  
  pdat1[["sd_value"]] <- pdat2[["value"]]
  pdat1[["variable"]] <- as.numeric(as.character(pdat1[["variable"]]))
  pdat1[["ci_min"]] <- pdat1[["value"]] - 2 * pdat1[["sd_value"]]
  pdat1[["ci_max"]] <- pdat1[["value"]] + 2 * pdat1[["sd_value"]]
  
  pdat1
}

pdat_1_1 <- process_result(abs(res_trajectory[[1]]), 
                           res_trajectory[[2]], 
                           row_name = u_method_type, 
                           col_name = trajectory_num_vec)
pdat_1_2 <- process_result(res_trajectory[[3]], 
                           res_trajectory[[4]], 
                           row_name = u_method_type, 
                           col_name = trajectory_num_vec)
pdat_2_1 <- process_result(abs(res_time[[1]]), res_time[[2]], 
                           row_name = u_method_type, 
                           col_name = trajectory_num_vec)
pdat_2_2 <- process_result(res_time[[3]], res_time[[4]], 
                           row_name = u_method_type, 
                           col_name = trajectory_num_vec)

pdat <- rbind.data.frame(cbind.data.frame(pdat_1_1, type = "logBias", class = "trajectory"), 
                         cbind.data.frame(pdat_1_2, type = "logMSE", class = "trajectory"), 
                         cbind.data.frame(pdat_2_1, type = "logBias", class = "time"), 
                         cbind.data.frame(pdat_2_2, type = "logMSE", class = "time"))

p <- ggplot(pdat, aes(x = variable, y = log10(value), color = Method)) + 
  facet_grid(type ~ class, scales = "free_y") + 
  geom_point(aes(group = Method), size = 2.5) +
  geom_line(aes(group = Method), size = 1) +
  geom_errorbar(aes(ymin = log10(ci_min), ymax = log10(ci_max)), size = 1, width = 1) + 
  scale_x_continuous(breaks = trajectory_num_vec) +
  theme_bw() + 
  xlab("") + ylab("") + 
  theme(legend.position = "bottom", 
        legend.box.margin = margin(t = -20, b = 10, r = 0, l = 0)) + 
  guides(guide_legend(byrow = TRUE))
p
ggsave(filename = "comparison.jpg", plot = p, width = 6, height = 5.6)

