###Last Function version 11/11/24


# x: Data frame containing the data to analyze.
# Replacement: Boolean (TRUE or FALSE). 
#               Indicates whether to sample with replacement (default = TRUE).
# GR: String specifying the grouping variable (e.g., "group"). 
#     This variable is used to define which observations belong to which groups.
# IN: Integer specifying the starting index for the traits/variables to analyze (default = 2).
# FIN: Integer specifying the ending index for the traits/variables to analyze (default = number of columns in x).
# G1: String specifying the first group to compare (default = "high").
# G2: String specifying the second group to compare (default = "low").
# n_iter: Integer specifying the number of bootstrap iterations (default = 1000).
# conf_method: String indicating the method to calculate confidence intervals. Options include:
#              "SD" for standard deviation, 
#              "percentile" for percentile method (default), 
#              "BCa" for bias-corrected and accelerated method, 
#              "bootstrap-t" for bootstrap-t method.
# SEED: Integer for setting the random seed for reproducibility (default = 123).
# method: String specifying the effect size calculation method. Options include:
#         "Hedges" for Hedges' g (default).
# show_histograms: Boolean (TRUE or FALSE) to indicate whether to display histograms of effect size distributions (default = TRUE).
# show_pairwise_comparisons: Can be:
#                            FALSE: No pairwise comparisons are shown.
#                            TRUE: Both original and absolute pairwise comparisons are shown.
#                            "ABS": Only absolute pairwise comparisons are shown.
#                            "PAIR": Only original pairwise comparisons are shown.

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
                        valid_pop1 <- pop1[[var]][!is.na(pop1[[var]])]
                        valid_pop2 <- pop2[[var]][!is.na(pop2[[var]])]
                        sample1 <- sample(valid_pop1, size = length(valid_pop1), replace = Replacement)
                        sample2 <- sample(valid_pop2, size = length(valid_pop2), replace = Replacement)
                        list(
                                effect = calc_effect_size(sample1, sample2, method),
                                n1 = length(sample1),
                                n2 = length(sample2)
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
}












results<-BootstrapEffectSizesOPTG(test_data,n_iter = 1000,method="Cohens", GR= "group",show_histograms = F,show_pairwise_comparisons=T)
results

#method="Cohens" and Hedges" should be clearly specified within the formula
##Pairwaise comparisons always performed. Shouyld be done only if show_pairwise_comparisons is NOT set to FALSE (or F). If set to false the pairwise comparison should not be done nor shown
##"ABS" and "PAIR" comparisons are not working as arguments of the formula. THey should be performed and shown if specified with this logic for the argument:
#T or TRUE: both normal and absolute values of the ES should pe compared and shown graphically and as a table in the results
#F or FALSE: none of the comparisons should be made nor shown 
#"ABS" only the absolute values should be compared and shown graphically and as a table in the results
#"PAIR" only the ormal values should be compared and shown graphically and as a table in the results
##Histograms should be not shown directly but stored in a object called Hist. It can be called if show_histograms=T. otherwise not shown.

##More profound changings to be done in the future:
#1. shorter arguments names in the function (e.g. "show_pairwise_comparisons" should be "PW_C"; "show_histograms" may be "S_H")
#2. Arguments G1 = "high" and  G2 = "low" are useless. The function should be reorganized in order to autonomously recognize the levels of the "GR=" arguments
#3. In order to make the funcion more general there should be a new argument. This argument is quite troublesome as it will change drastically the function and may need a complete new chuck of text:
#the argument should be called "RPS=" (meaning "Resampled Population Size"). IF present it may have 2 kind of values: 
# TRUE (or T) the function will generate the new populations (as always excluding the NA when resampling and when calculating the new population size) based on the smallest pupulation size based on the traits presents in the dataset (i.e. the population size is based for all the traits on the smallest among the traits)
# a number ranging from 2 to N. In this case for all the traits the specified number will be used as a resampled populations size. AS always the function should be able to deal the NAs always excluding them (there should be no NA in the resampled datasets).
# Expedially for this last argument, in the "Parameters used:" section the used arguments and meaning of the function has to be clear.