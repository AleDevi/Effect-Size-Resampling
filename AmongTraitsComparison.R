#######Updated version to compare effect sizes among traits visually
BootstrapEffectSizesOPTG <- function(x, Replacement = TRUE, GR = "group", IN = 2, FIN = ncol(x), G1 = "high", G2 = "low", 
                                     n_iter = 1000, conf_method = "percentile", SEED = 123, method = "Hedges") {
        set.seed(SEED)
        library(dplyr)
        library(ggplot2)
        library(tidyr)
        library(boot)  # Necessario per il metodo BCa
        
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
        
        effect_sizes <- replicate(n_iter, {
                sapply(vars, function(var) {
                        sample1 <- sample(pop1[[var]], size = sum(!is.na(pop1[[var]])), replace = Replacement)
                        sample2 <- sample(pop2[[var]], size = sum(!is.na(pop2[[var]])), replace = Replacement)
                        calc_effect_size(sample1, sample2, method)
                })
        })
        
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
                               ci <- try(boot.ci(boot_obj, type = "stud")$student[4:5], silent = TRUE)
                               if (inherits(ci, "try-error")) quantile(x, c(0.025, 0.975), na.rm = TRUE) else ci
                       })
        }
        
        ci_summary <- data.frame(
                Variable = vars,
                mean = rowMeans(effect_sizes, na.rm = TRUE)
        )
        
        ci_results <- t(sapply(1:nrow(effect_sizes), function(i) calc_ci(effect_sizes[i,], conf_method)))
        ci_summary$lower <- ci_results[,1]
        ci_summary$upper <- ci_results[,2]
        
        ci_summary <- ci_summary %>%
                mutate(
                        significant = lower > 0 | upper < 0,
                        observed_n1 = sapply(vars, function(var) sum(!is.na(pop1[[var]]))),
                        observed_n2 = sapply(vars, function(var) sum(!is.na(pop2[[var]]))),
                        mean_resampled_n1 = sapply(vars, function(var) sum(!is.na(pop1[[var]]))),
                        mean_resampled_n2 = sapply(vars, function(var) sum(!is.na(pop2[[var]])))
                )
        
        # Plotting the effect sizes
        plot <- ggplot(ci_summary, aes(x = Variable, y = mean, ymin = lower, ymax = upper)) +
                geom_pointrange(color = "blue", size = 1) +
                geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
                labs(title = "Effect Sizes with Confidence Intervals", x = "Variables", y = "Effect Size") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        print(plot)
        cat("\nMethod used for Effect Sizes:", method, "\n")
        cat("Confidence Method:", conf_method, "\n")
        cat("Number of Iterations:", n_iter, "\n\n")
        
        histograms <- lapply(1:nrow(effect_sizes), function(i) {
                ggplot(data.frame(effect_dist = effect_sizes[i,]), aes(x = effect_dist)) +
                        geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
                        geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
                        labs(title = paste("Effect Size Distribution for", vars[i]), x = "Effect Size", y = "Frequency") +
                        theme_minimal()
        })
        
        # Pairwise comparisons among traits
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
        
        # Plot for pairwise comparisons
        pairwise_plot <- ggplot(pairwise_comparisons, aes(x = Trait1, y = Trait2, fill = Significant_Difference)) +
                geom_tile(color = "white") +
                scale_fill_manual(values = c("TRUE" = "green", "FALSE" = "grey")) +
                labs(title = "Pairwise Significant Differences Between Traits", x = "Trait 1", y = "Trait 2") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        print(pairwise_plot)
        
        return(list(summary = ci_summary, pairwise_comparisons = pairwise_comparisons, histograms = histograms))
}
