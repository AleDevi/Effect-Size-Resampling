###BootstrapEffectSizes (BES). a new version obtained and tuned with blackboxAi
#this version should have all the features of those before 

BootstrapEffectSizesOPTG <- function(x, Replacement = TRUE, GR = "group", IN = 2, FIN = ncol(x), G1 = "high", G2 = "low", 
                                     n_iter = 1000, conf_method = "percentile", SEED = 123, method = "Hedges",
                                     show_histograms = TRUE, show_pairwise_comparisons = FALSE) {
        
        set.seed(SEED)
        library(dplyr)
        library(ggplot2)
        library(tidyr)
        library(boot)  # Necessary for the BCa method
        
        calc_effect_size <- function(pop1, pop2, method) {
                pop1 <- pop1[!is.na(pop1)]
                pop2 <- pop2[!is.na(pop2)]
                
                n1 <- length(pop1)
                n2 <- length(pop2)
                
                if (n1 < 2 || n2 < 2) return(NA)
                
                p1mean <- mean(pop1)
                p2mean <- mean(pop2)
                p1sd <- sd(pop1)
                p2sd <- sd(pop2)
                
                pooled_sd <- sqrt(((n1 - 1) * p1sd^2 + (n2 - 1) * p2sd^2) / (n1 + n2 - 2))
                effect <- (p1mean - p2mean) / pooled_sd
                
                if (method == "Hedges") {
                        effect <- effect * (1 - (3 / (4 * (n1 + n2) - 9)))
                }
                
                return(effect)
        }
        
        vars <- names(x)[IN:FIN]
        pop1 <- x[x[[GR]] == G1, vars, drop = FALSE]
        pop2 <- x[x[[GR]] == G2, vars, drop = FALSE]
        
        # Calcolo degli effect sizes e tracciamento del numero di campioni
        effect_sizes <- replicate(n_iter, {
                sapply(vars, function(var) {
                        sample1 <- sample(pop1[[var]], size = sum(!is.na(pop1[[var]])), replace = Replacement)
                        sample2 <- sample(pop2[[var]], size = sum(!is.na(pop2[[var]])), replace = Replacement)
                        list(
                                effect = calc_effect_size(sample1, sample2, method),
                                n1 = sum(!is.na(sample1)),
                                n2 = sum(!is.na(sample2))
                        )
                }, simplify = FALSE)
        }, simplify = FALSE)
        
        # Estrazione degli effect sizes e dei numeri di campioni
        effect_sizes_matrix <- sapply(effect_sizes, function(iter) sapply(iter, function(var) var$effect))
        n1_matrix <- sapply(effect_sizes, function(iter) sapply(iter, function(var) var$n1))
        n2_matrix <- sapply(effect_sizes, function(iter) sapply(iter, function(var) var$n2))
        
        calc_ci <- function(x, conf_method) {
                switch(conf_method,
                       "SD" = c(mean(x, na.rm = TRUE) - sd(x, na.rm = TRUE),
                                mean(x, na.rm = TRUE) + sd(x, na.rm = TRUE)),
                       "percentile" = quantile(x, c(0.025, 0.975), na.rm = TRUE),
                       "BCa" = {
                               boot_obj <- boot(x, function(d, i) mean(d[i]), R = 1000)
                               ci <- try(boot.ci(boot_obj, type = "bca")$bca[4:5], silent = TRUE)
                               if (inherits(ci, "try-error")) quantile(x, c(0.025, 0.975), na.rm = TRUE) else ci
                       },
                       "bootstrap-t" = {
                               boot_obj <- boot(x, function(d, i) c(mean(d[i]), sd(d[i]) / sqrt(length(d[i]))), R = 1000)
                               ci <- try(boot.ci( boot_obj, type = "stud")$student[4:5], silent = TRUE)
                               if (inherits(ci, "try-error")) quantile(x, c(0.025, 0.975), na.rm = TRUE) else ci
                       })
        }
        
        # Calculate means and confidence intervals for each trait
        ci_summary <- data.frame(
                Variable = vars,
                mean = rowMeans(effect_sizes_matrix, na.rm = TRUE)
        )
        
        ci_results <- t(sapply(1:nrow(effect_sizes_matrix), function(i) calc_ci(effect_sizes_matrix[i,], conf_method)))
        ci_summary$lower <- ci_results[, 1]
        ci_summary$upper <- ci_results[, 2]
        
        # Add observed sample sizes and mean resampling numbers
        ci_summary$observed_n1 <- sapply(vars, function(var) sum(!is.na(pop1[[var]])))
        ci_summary$observed_n2 <- sapply(vars, function(var) sum(!is.na(pop2[[var]])))
        ci_summary$mean_resampled_n1 <- rowMeans(n1_matrix)
        ci_summary$mean_resampled_n2 <- rowMeans(n2_matrix)
        
        # Plot effect sizes with confidence intervals
        plot <- ggplot(ci_summary, aes(x = Variable, y = mean, ymin = lower, ymax = upper)) +
                geom_pointrange(color = "blue", size = 1) +
                geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
                labs(title = "Effect Sizes with Confidence Intervals", x = "Variables", y = "Effect Size") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        print(plot)
        
        # Calculate absolute effect sizes and confidence intervals
        abs_effect_sizes <- abs(effect_sizes_matrix)
        abs_ci_summary <- data.frame(
                Variable = vars,
                abs_mean = rowMeans(abs_effect_sizes, na.rm = TRUE)
        )
        
        abs_ci_results <- t(sapply(1:nrow(abs_effect_sizes), function(i) calc_ci(abs_effect_sizes[i,], conf_method)))
        abs_ci_summary <- abs_ci_summary %>%
                mutate(lower_abs = abs_ci_results[, 1], upper_abs = abs_ci_results[, 2])
        
        # Generate pairwise comparisons
        pairwise_comparisons <- expand.grid(Trait1 = vars, Trait2 = vars) %>%
                filter(Trait1 != Trait2) %>%
                rowwise() %>%
                mutate(
                        lower1 = ci_summary$lower[ci_summary$Variable == Trait1],
                        upper1 = ci_summary$upper[ci_summary$Variable == Trait1],
                        lower2 = ci_summary$lower[ci_summary$Variable == Trait2],
                        upper2 = ci_summary$upper[ci_summary$Variable == Trait2],
                        Significant_Difference = !(upper1 > lower2 & lower1 < upper2)
                ) %>%
                ungroup() %>%
                select(Trait1, Trait2, Significant_Difference)
        
        # Absolute pairwise comparisons
        abs_pairwise_comparisons <- expand.grid(Trait1 = vars, Trait2 = vars) %>%
                filter(Trait1 != Trait2) %>%
                rowwise() %>%
                mutate(
                        lower1_abs = abs_ci_summary$lower_abs[abs_ci_summary$Variable == Trait1],
                        upper1_abs = abs_ci_summary$upper_abs[abs_ci_summary$Variable == Trait1],
                        lower2_abs = abs_ci_summary$lower_abs[abs_ci_summary$Variable == Trait2],
                        upper2_abs = abs_ci_summary$upper_abs[abs_ci_summary$Variable == Trait2],
                        Abs_Significant_Difference = !(upper1_abs > lower2_abs & lower1_abs < upper2_abs)
                ) %>%
                ungroup() %>%
                select(Trait1, Trait2, Abs_Significant_Difference)
        
        # Show histograms if requested
        if (show_histograms) {
                histograms <- lapply(1:nrow(effect_sizes_matrix), function(i) {
                        ggplot(data.frame(effect_dist = effect_sizes_matrix[i,]), aes(x = effect_dist)) +
                                geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
                                geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
                                labs(title = paste("Effect Size Distribution for", vars[i]), x = "Effect Size", y = "Frequency") +
                                theme_minimal()
                })
                
                # Print histograms
                for (hist in histograms) {
                        print(hist)
                }
        }
        
        # Show pairwise comparisons based on the argument
        if (show_pairwise_comparisons) {
                # Visualization of pairwise comparisons on original values
                pairwise_plot <- ggplot(pairwise_comparisons, aes(x = Trait1, y = Trait2, fill = Significant_Difference)) +
                        geom_tile(color = "white") +
                        scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "grey")) +
                        labs(title = "Pairwise Significant Differences Between Traits (Original ES)", x = "Trait 1", y = "Trait 2") +
                        theme_minimal() +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                print(pairwise_plot)
                
                # Visualization of pairwise comparisons on absolute values
                abs_pairwise_plot <- ggplot(abs_pairwise_comparisons, aes(x = Trait1, y = Trait2, fill = Abs_Significant_Difference)) +
                        geom_tile(color = "white") +
                        scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "grey")) +
                        labs(title = "Pairwise Significant Differences Between Traits (Absolute ES)", x = "Trait 1", y = "Trait 2") +
                        theme_minimal() +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                print(abs_pairwise_plot)
        }
        
        # Print parameter details
        cat("\nParameters used:\n")
        cat("Number of permutations (n_iter):", n_iter, "\n")
        cat("Type of effect size calculated (method):", method, "\n")
        cat("Confidence method used:", conf_method, "\n")
        cat("Group variable (GR):", GR, "\n")
        cat("Group 1 (G1):", G1, "\n")
        cat("Group 2 (G2):", G2, "\n")
        cat("Replacement used in sampling:", Replacement, "\n")
        
        return(list(summary = ci_summary, pairwise_comparisons = pairwise_comparisons, abs_pairwise_comparisons = abs_pairwise_comparisons))
}#sbaglia l'andling degli NA

#########################################################
############################################################
#########################################################
############################################################
#########################################################
############################################################


BootstrapEffectSizesOPTG <- function(x, Replacement = TRUE, GR = "group", IN = 2, FIN = ncol(x), G1 = "high", G2 = "low", 
                                     n_iter = 1000, conf_method = "percentile", SEED = 123, method = "Hedges",
                                     show_histograms = TRUE, show_pairwise_comparisons = FALSE) {
        
        set.seed(SEED)
        library(dplyr)
        library(ggplot2)
        library(tidyr)
        library(boot)  # Necessary for the BCa method
        
        # Effect size calculation function
        calc_effect_size <- function(pop1, pop2, method) {
                if (length(pop1) < 2 || length(pop2) < 2) return(NA)
                
                p1mean <- mean(pop1)
                p2mean <- mean(pop2)
                pooled_sd <- sqrt(((length(pop1) - 1) * var(pop1) + (length(pop2) - 1) * var(pop2)) / (length(pop1) + length(pop2) - 2))
                effect <- (p1mean - p2mean) / pooled_sd
                
                if (method == "Hedges") {
                        effect <- effect * (1 - (3 / (4 * (length(pop1) + length(pop2)) - 9)))
                }
                
                return(effect)
        }
        
        # Extract variable names for analysis
        vars <- names(x)[IN:FIN]
        
        # Prepare non-NA samples for each trait independently
        pop1 <- lapply(vars, function(var) x[[var]][x[[GR]] == G1 & !is.na(x[[var]])])
        pop2 <- lapply(vars, function(var) x[[var]][x[[GR]] == G2 & !is.na(x[[var]])])
        names(pop1) <- names(pop2) <- vars
        
        # Track observed sample sizes (non-NA counts) for each trait
        observed_n1 <- sapply(pop1, length)
        observed_n2 <- sapply(pop2, length)
        
        # Perform bootstrap resampling and calculate effect sizes
        effect_sizes <- replicate(n_iter, {
                sapply(vars, function(var) {
                        sample1 <- sample(pop1[[var]], size = observed_n1[var], replace = Replacement)
                        sample2 <- sample(pop2[[var]], size = observed_n2[var], replace = Replacement)
                        calc_effect_size(sample1, sample2, method)
                })
        })
        
        # Summarize results
        ci_summary <- data.frame(
                Variable = vars,
                mean = rowMeans(effect_sizes, na.rm = TRUE)
        )
        
        # Calculate confidence intervals
        ci_results <- t(apply(effect_sizes, 1, function(es) {
                switch(conf_method,
                       "SD" = c(mean(es, na.rm = TRUE) - sd(es, na.rm = TRUE), mean(es, na.rm = TRUE) + sd(es, na.rm = TRUE)),
                       "percentile" = quantile(es, c(0.025, 0.975), na.rm = TRUE),
                       "BCa" = {
                               boot_obj <- boot(es, function(d, i) mean(d[i]), R = 1000)
                               ci <- try(boot.ci(boot_obj, type = "bca")$bca[4:5], silent = TRUE)
                               if (inherits(ci, "try-error")) quantile(es, c(0.025, 0.975), na.rm = TRUE) else ci
                       },
                       "bootstrap-t" = {
                               boot_obj <- boot(es, function(d, i) c(mean(d[i]), sd(d[i]) / sqrt(length(d[i]))), R = 1000)
                               ci <- try(boot.ci(boot_obj, type = "stud")$student[4:5], silent = TRUE)
                               if (inherits(ci, "try-error")) quantile(es, c(0.025, 0.975), na.rm = TRUE) else ci
                       })
        }))
        
        ci_summary$lower <- ci_results[, 1]
        ci_summary$upper <- ci_results[, 2]
        ci_summary$observed_n1 <- observed_n1
        ci_summary$observed_n2 <- observed_n2
        ci_summary$mean_resampled_n1 <- observed_n1
        ci_summary$mean_resampled_n2 <- observed_n2
        
        # Plot effect sizes with confidence intervals
        plot <- ggplot(ci_summary, aes(x = Variable, y = mean, ymin = lower, ymax = upper)) +
                geom_pointrange(color = "blue", size = 1) +
                geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
                labs(title = "Effect Sizes with Confidence Intervals", x = "Variables", y = "Effect Size") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        print(plot)
        
        # Show histograms if requested
        if (show_histograms) {
                histograms <- lapply(1:nrow(effect_sizes), function(i) {
                        ggplot(data.frame(effect_dist = effect_sizes[i,]), aes(x = effect_dist)) +
                                geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
                                geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
                                labs(title = paste("Effect Size Distribution for", vars[i]), x = "Effect Size", y = "Frequency") +
                                theme_minimal()
                })
                
                for (hist in histograms) {
                        print(hist)
                }
        }
        
        return(list(summary = ci_summary))
}###mancano le comparazioni pairwise
#########################################################
############################################################
#########################################################
############################################################
#########################################################
############################################################

BootstrapEffectSizesOPTG <- function(x, Replacement = TRUE, GR = "group", IN = 2, FIN = ncol(x), G1 = "high", G2 = "low", 
                                     n_iter = 1000, conf_method = "percentile", SEED = 123, method = "Hedges",
                                     show_histograms = TRUE, show_pairwise_comparisons = FALSE) {
        
        set.seed(SEED)
        library(dplyr)
        library(ggplot2)
        library(tidyr)
        library(boot)  # Necessary for the BCa method
        
        # Effect size calculation function
        calc_effect_size <- function(pop1, pop2, method) {
                if (length(pop1) < 2 || length(pop2) < 2) return(NA)
                
                p1mean <- mean(pop1)
                p2mean <- mean(pop2)
                pooled_sd <- sqrt(((length(pop1) - 1) * var(pop1) + (length(pop2) - 1) * var(pop2)) / (length(pop1) + length(pop2) - 2))
                effect <- (p1mean - p2mean) / pooled_sd
                
                if (method == "Hedges") {
                        effect <- effect * (1 - (3 / (4 * (length(pop1) + length(pop2)) - 9)))
                }
                
                return(effect)
        }
        
        # Extract variable names for analysis
        vars <- names(x)[IN:FIN]
        
        # Prepare non-NA samples for each trait independently
        pop1 <- lapply(vars, function(var) x[[var]][x[[GR]] == G1 & !is.na(x[[var]])])
        pop2 <- lapply(vars, function(var) x[[var]][x[[GR]] == G2 & !is.na(x[[var]])])
        names(pop1) <- names(pop2) <- vars
        
        # Track observed sample sizes (non-NA counts) for each trait
        observed_n1 <- sapply(pop1, length)
        observed_n2 <- sapply(pop2, length)
        
        # Perform bootstrap resampling and calculate effect sizes
        effect_sizes <- replicate(n_iter, {
                sapply(vars, function(var) {
                        sample1 <- sample(pop1[[var]], size = observed_n1[var], replace = Replacement)
                        sample2 <- sample(pop2[[var]], size = observed_n2[var], replace = Replacement)
                        calc_effect_size(sample1, sample2, method)
                })
        })
        
        # Summarize results
        ci_summary <- data.frame(
                Variable = vars,
                mean = rowMeans(effect_sizes, na.rm = TRUE)
        )
        
        # Calculate confidence intervals
        ci_results <- t(apply(effect_sizes, 1, function(es) {
                switch(conf_method,
                       "SD" = c(mean(es, na.rm = TRUE) - sd(es, na.rm = TRUE), mean(es, na.rm = TRUE) + sd(es, na.rm = TRUE)),
                       "percentile" = quantile(es, c(0.025, 0.975), na.rm = TRUE),
                       "BCa" = {
                               boot_obj <- boot(es, function(d, i) mean(d[i]), R = 1000)
                               ci <- try(boot.ci(boot_obj, type = "bca")$bca[4:5], silent = TRUE)
                               if (inherits(ci, "try-error")) quantile(es, c(0.025, 0.975), na.rm = TRUE) else ci
                       },
                       "bootstrap-t" = {
                               boot_obj <- boot(es, function(d, i) c(mean(d[i]), sd(d[i]) / sqrt(length(d[i]))), R = 1000)
                               ci <- try(boot.ci(boot_obj, type = "stud")$student[4:5], silent = TRUE)
                               if (inherits(ci, "try-error")) quantile(es, c(0.025, 0.975), na.rm = TRUE) else ci
                       })
        }))
        
        ci_summary$lower <- ci_results[, 1]
        ci_summary$upper <- ci_results[, 2]
        ci_summary$observed_n1 <- observed_n1
        ci_summary$observed_n2 <- observed_n2
        ci_summary$mean_resampled_n1 <- observed_n1
        ci_summary$mean_resampled_n2 <- observed_n2
        
        # Plot effect sizes with confidence intervals
        plot <- ggplot(ci_summary, aes(x = Variable, y = mean, ymin = lower, ymax = upper)) +
                geom_pointrange(color = "blue", size = 1) +
                geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
                labs(title = "Effect Sizes with Confidence Intervals", x = "Variables", y = "Effect Size") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        print(plot)
        
        # Calculate absolute effect sizes and confidence intervals
        abs_effect_sizes <- abs(effect_sizes)
        abs_ci_summary <- data.frame(
                Variable = vars,
                abs_mean = rowMeans(abs_effect_sizes, na.rm = TRUE)
        )
        
        abs_ci_results <- t(apply(abs_effect_sizes, 1, function(es) {
                switch(conf_method,
                       "SD" = c(mean(es, na.rm = TRUE) - sd(es, na.rm = TRUE), mean(es, na.rm = TRUE) + sd(es, na.rm = TRUE)),
                       "percentile" = quantile(es, c(0.025, 0.975), na.rm = TRUE),
                       "BCa" = {
                               boot_obj <- boot(es, function(d, i) mean(d[i]), R = 1000)
                               ci <- try(boot.ci(boot_obj, type = "bca")$bca[4:5], silent = TRUE)
                               if (inherits(ci, "try-error")) quantile(es, c(0.025, 0.975), na.rm = TRUE) else ci
                       },
                       "bootstrap-t" = {
                               boot_obj <- boot(es, function(d, i) c(mean(d[i]), sd(d[i]) / sqrt(length(d[i]))), R = 1000)
                               ci <- try(boot.ci(boot_obj, type = "stud")$student[4:5], silent = TRUE)
                               if (inherits(ci, "try-error")) quantile(es, c(0.025, 0.975), na.rm = TRUE) else ci
                       })
        }))
        
        abs_ci_summary <- abs_ci_summary %>%
                mutate(lower_abs = abs_ci_results[, 1], upper_abs = abs_ci_results[, 2])
        
        # Generate pairwise comparisons
        pairwise_comparisons <- expand.grid(Trait1 = vars, Trait2 = vars) %>%
                filter(Trait1 != Trait2) %>%
                rowwise() %>%
                mutate(
                        lower1 = ci_summary$lower[ci_summary$Variable == Trait1],
                        upper1 = ci_summary$upper[ci_summary$Variable == Trait1],
                        lower2 = ci_summary$lower[ci_summary$Variable == Trait2],
                        upper2 = ci_summary$upper[ci_summary$Variable == Trait2],
                        Significant_Difference = !(upper1 > lower2 & lower1 < upper2)
                ) %>%
                ungroup() %>%
                select(Trait1, Trait2, Significant_Difference)
        
        # Absolute pairwise comparisons
        abs_pairwise_comparisons <- expand.grid(Trait1 = vars, Trait2 = vars) %>%
                filter(Trait1 != Trait2) %>%
                rowwise() %>%
                mutate(
                        lower1_abs = abs_ci_summary$lower_abs[abs_ci_summary$Variable == Trait1],
                        upper1_abs = abs_ci_summary$upper_abs[abs_ci_summary$Variable == Trait1],
                        lower2_abs = abs_ci_summary$lower_abs[abs_ci_summary$Variable == Trait2],
                        upper2_abs = abs_ci_summary$upper_abs[abs_ci_summary$Variable == Trait2],
                        Abs_Significant_Difference = !(upper1_abs > lower2_abs & lower1_abs < upper2_abs)
                ) %>%
                ungroup() %>%
                select(Trait1, Trait2, Abs_Significant_Difference)
        
        # Show histograms if requested
        if (show_histograms) {
                histograms <- lapply(1:nrow(effect_sizes), function(i) {
                        ggplot(data.frame(effect_dist = effect_sizes[i,]), aes(x = effect_dist)) +
                                geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
                                geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
                                labs(title = paste("Effect Size Distribution for", vars[i]), x = "Effect Size", y = "Frequency") +
                                theme_minimal()
                })
                
                for (hist in histograms) {
                        print(hist)
                }
        }
        
        # Show pairwise comparisons if requested
        if (show_pairwise_comparisons == TRUE) {
                print(pairwise_comparisons)
        } else if (show_pairwise_comparisons == "ABS") {
                print(abs_pairwise_comparisons)
        } else if (show_pairwise_comparisons == "PAIR") {
                print(pairwise_comparisons)
                print(abs_pairwise_comparisons)
        }
        
        return(list(summary = ci_summary, pairwise_comparisons = pairwise_comparisons, abs_pairwise_comparisons = abs_pairwise_comparisons))
}###sembra OK ma manca la parte grafica

################################################################
#################################################################
##############################################################
###############################################################

BootstrapEffectSizesOPTG <- function(x, Replacement = TRUE, GR = "group", IN = 2, FIN = ncol(x), G1 = "high", G2 = "low", 
                                     n_iter = 1000, conf_method = "percentile", SEED = 123, method = "Hedges",
                                     show_histograms = TRUE, show_pairwise_comparisons = FALSE) {
        
        set.seed(SEED)
        library(dplyr)
        library(ggplot2)
        library(tidyr)
        library(boot)  # Necessary for the BCa method
        
        # Effect size calculation function
        calc_effect_size <- function(pop1, pop2, method) {
                if (length(pop1) < 2 || length(pop2) < 2) return(NA)
                
                p1mean <- mean(pop1)
                p2mean <- mean(pop2)
                pooled_sd <- sqrt(((length(pop1) - 1) * var(pop1) + (length(pop2) - 1) * var(pop2)) / (length(pop1) + length(pop2) - 2))
                effect <- (p1mean - p2mean) / pooled_sd
                
                if (method == "Hedges") {
                        effect <- effect * (1 - (3 / (4 * (length(pop1) + length(pop2)) - 9)))
                }
                
                return(effect)
        }
        
        # Extract variable names for analysis
        vars <- names(x)[IN:FIN]
        
        # Prepare non-NA samples for each trait independently
        pop1 <- lapply(vars, function(var) x[[var]][x[[GR]] == G1 & !is.na(x[[var]])])
        pop2 <- lapply(vars, function(var) x[[var]][x[[GR]] == G2 & !is.na(x[[var]])])
        names(pop1) <- names(pop2) <- vars
        
        # Track observed sample sizes (non-NA counts) for each trait
        observed_n1 <- sapply(pop1, length)
        observed_n2 <- sapply(pop2, length)
        
        # Perform bootstrap resampling and calculate effect sizes
        effect_sizes <- replicate(n_iter, {
                sapply(vars, function(var) {
                        sample1 <- sample(pop1[[var]], size = observed_n1[var], replace = Replacement)
                        sample2 <- sample(pop2[[var]], size = observed_n2[var], replace = Replacement)
                        calc_effect_size(sample1, sample2, method)
                })
        })
        
        # Summarize results
        ci_summary <- data.frame(
                Variable = vars,
                mean = rowMeans(effect_sizes, na.rm = TRUE)
        )
        
        # Calculate confidence intervals
        ci_results <- t(apply(effect_sizes, 1, function(es) {
                switch(conf_method,
                       "SD" = c(mean(es, na.rm = TRUE) - sd(es, na.rm = TRUE), mean(es, na.rm = TRUE) + sd(es, na.rm = TRUE)),
                       "percentile" = quantile(es, c(0.025, 0.975), na.rm = TRUE),
                       "BCa" = {
                               boot_obj <- boot(es, function(d, i) mean(d[i]), R = 1000)
                               ci <- try(boot.ci(boot_obj, type = "bca")$bca[4:5], silent = TRUE)
                               if (inherits(ci, "try-error")) quantile(es, c(0.025, 0.975), na.rm = TRUE) else ci
                       },
                       "bootstrap-t" = {
                               boot_obj <- boot(es, function(d, i) c(mean(d[i]), sd(d[i]) / sqrt(length(d[i]))), R = 1000)
                               ci <- try(boot.ci(boot_obj, type = "stud")$student[4:5], silent = TRUE)
                               if (inherits(ci, "try-error")) quantile(es, c(0.025, 0.975), na.rm = TRUE) else ci
                       })
        }))
        
        ci_summary$lower <- ci_results[, 1]
        ci_summary$upper <- ci_results[, 2]
        ci_summary$observed_n1 <- observed_n1
        ci_summary$observed_n2 <- observed_n2
        ci_summary$mean_resampled_n1 <- observed_n1
        ci_summary$mean_resampled_n2 <- observed_n2
        
        # Plot effect sizes with confidence intervals
        plot <- ggplot(ci_summary, aes(x = Variable, y = mean, ymin = lower, ymax = upper)) +
                geom_pointrange(color = "blue", size = 1) +
                geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
                labs(title = "Effect Sizes with Confidence Intervals", x = "Variables", y = "Effect Size") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        print(plot)
        
        # Calculate absolute effect sizes and confidence intervals
        abs_effect_sizes <- abs(effect_sizes)
        abs_ci_summary <- data.frame(
                Variable = vars,
                abs_mean = rowMeans(abs_effect_sizes, na.rm = TRUE)
        )
        
        abs_ci_results <- t(apply(abs_effect_sizes, 1, function(es) {
                switch(conf_method,
                       "SD" = c(mean(es, na.rm = TRUE) - sd(es, na.rm = TRUE), mean(es, na.rm = TRUE) + sd(es, na.rm = TRUE)),
                       "percentile" = quantile(es, c(0.025, 0.975), na.rm = TRUE),
                       "BCa" = {
                               boot_obj <- boot(es, function(d, i) mean(d[i]), R = 1000)
                               ci <- try(boot.ci(boot_obj, type = "bca")$bca[4:5], silent = TRUE)
                               if (inherits(ci, "try-error")) quantile(es, c(0.025, 0.975), na.rm = TRUE) else ci
                       },
                       "bootstrap-t" = {
                               boot_obj <- boot(es, function(d, i) c(mean(d[i]), sd(d[i]) / sqrt(length(d[i]))), R = 1000)
                               ci <- try(boot.ci(boot_obj, type = "stud")$student[4:5], silent = TRUE)
                               if (inherits(ci, "try-error")) quantile(es, c(0.025, 0.975), na.rm = TRUE) else ci
                       })
        }))
        
        abs_ci_summary <- abs_ci_summary %>%
                mutate(lower_abs = abs_ci_results[, 1], upper_abs = abs_ci_results[, 2])
        
        # Generate pairwise comparisons
        pairwise_comparisons <- expand.grid(Trait1 = vars, Trait2 = vars) %>%
                filter(Trait1 != Trait2) %>%
                rowwise() %>%
                mutate(
                        lower1 = ci_summary$lower[ci_summary$Variable == Trait1],
                        upper1 = ci_summary$upper[ci_summary$Variable == Trait1],
                        lower2 = ci_summary$lower[ci_summary$Variable == Trait2],
                        upper2 = ci_summary$upper[ci_summary$Variable == Trait2],
                        Significant_Difference = !(upper1 > lower2 & lower1 < upper2)
                ) %>%
                ungroup() %>%
                select(Trait1, Trait2, Significant_Difference)
        
        # Absolute pairwise comparisons
        abs_pairwise_comparisons <- expand.grid(Trait1 = vars, Trait2 = vars) %>%
                filter(Trait1 != Trait2) %>%
                rowwise() %>%
                mutate(
                        lower1_abs = abs_ci_summary$lower_abs[abs_ci_summary$Variable == Trait1],
                        upper1_abs = abs_ci_summary$upper_abs[abs_ci_summary$Variable == Trait1],
                        lower2_abs = abs_ci_summary$lower_abs[abs_ci_summary$Variable == Trait2],
                        upper2_abs = abs_ci_summary$upper_abs[abs_ci_summary$Variable == Trait2],
                        Abs_Significant_Difference = !(upper1_abs > lower2_abs & lower1_abs < upper2_abs)
                ) %>%
                ungroup() %>%
                select(Trait1, Trait2, Abs_Significant_Difference)
        
        # Show histograms if requested
        if (show_histograms) {
                histograms <- lapply(1:nrow(effect_sizes), function(i) {
                        ggplot(data.frame(effect_dist = effect_sizes[i,]), aes(x = effect_dist)) +
                                geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
                                geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
                                labs(title = paste("Effect Size Distribution for", vars[i]), x = "Effect Size", y = "Frequency") +
                                theme_minimal()
                })
                
                for (hist in histograms) {
                        print(hist)
                }
        }
        
        # Show pairwise comparisons in heatmaps based on input argument
        if (show_pairwise_comparisons == TRUE) {
                # Pairwise heatmap for original effect sizes
                pairwise_plot <- ggplot(pairwise_comparisons, aes(x = Trait1, y = Trait2, fill = Significant_Difference)) +
                        geom_tile(color = "white") +
                        scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "grey")) +
                        labs(title = "Pairwise Significant Differences Between Traits (Original ES)", x = "Trait 1", y = "Trait 2") +
                        theme_minimal() +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                print(pairwise_plot)
                
                # Pairwise heatmap for absolute effect sizes
                abs_pairwise_plot <- ggplot(abs_pairwise_comparisons, aes(x = Trait1, y = Trait2, fill = Abs_Significant_Difference)) +
                        geom_tile(color = "white") +
                        scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "grey")) +
                        labs(title = "Pairwise Significant Differences Between Traits (Absolute ES)", x = "Trait 1", y = "Trait 2") +
                        theme_minimal() +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                print(abs_pairwise_plot)
                
        } else if (show_pairwise_comparisons == "PAIR") {
                # Pairwise heatmap for original effect sizes
                pairwise_plot <- ggplot(pairwise_comparisons, aes(x = Trait1, y = Trait2, fill = Significant_Difference)) +
                        geom_tile(color = "white") +
                        scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "grey")) +
                        labs(title = "Pairwise Significant Differences Between Traits (Original ES)", x = "Trait 1", y = "Trait 2") +
                        theme_minimal() +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                print(pairwise_plot)
                
        } else if (show_pairwise_comparisons == "ABS") {
                # Pairwise heatmap for absolute effect sizes
                abs_pairwise_plot <- ggplot(abs_pairwise_comparisons, aes(x = Trait1, y = Trait2, fill = Abs_Significant_Difference)) +
                        geom_tile(color = "white") +
                        scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "grey")) +
                        labs(title = "Pairwise Significant Differences Between Traits (Absolute ES)", x = "Trait 1", y = "Trait 2") +
                        theme_minimal() +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                print(abs_pairwise_plot)
        }
        
        return(list(summary = ci_summary, pairwise_comparisons = pairwise_comparisons, abs_pairwise_comparisons = abs_pairwise_comparisons))
}






results<-BootstrapEffectSizesOPTG(test_data,n_iter = 100,method="Cohens", GR= "group",show_histograms = F, show_pairwise_comparisons = F)
results
