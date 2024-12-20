###BootstrapEffectSizesOPTG the last verision of the function


BootstrapEffectSizesOPTG <- function(x, Replacement = TRUE, GR = "group", IN = 2, FIN = ncol(x), G1 = "high", G2 = "low", 
                                    n_iter = 1000, conf_method = "percentile", SEED = 123, method = "Hedges",
                                    show_histograms = TRUE, show_pairwise_comparisons = FALSE) {
       
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
        
        # Calcola le medie e gli intervalli di confidenza originali per ciascun tratto
        ci_summary <- data.frame(
                Variable = vars,
                mean = rowMeans(effect_sizes, na.rm = TRUE)
        )
        
        ci_results <- t(sapply(1:nrow(effect_sizes), function(i) calc_ci(effect_sizes[i,], conf_method)))
        ci_summary$lower <- ci_results[, 1]
        ci_summary$upper <- ci_results[, 2]
        
        # Plot degli effect sizes con gli intervalli di confidenza (senza valore assoluto)
        plot <- ggplot(ci_summary, aes(x = Variable, y = mean, ymin = lower, ymax = upper)) +
                geom_pointrange(color = "blue", size = 1) +
                geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
                labs(title = "Effect Sizes with Confidence Intervals", x = "Variables", y = "Effect Size") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        print(plot)
        
        # Calcola gli intervalli di confidenza per i valori assoluti per il confronto tra ES
        abs_effect_sizes <- abs(effect_sizes)
        abs_ci_summary <- data.frame(
                Variable = vars,
                abs_mean = rowMeans(abs_effect_sizes, na.rm = TRUE)
        )
        
        abs_ci_results <- t(sapply(1:nrow(abs_effect_sizes), function(i) calc_ci(abs_effect_sizes[i,], conf_method)))
        abs_ci_summary <- abs_ci_summary %>%
                mutate(lower_abs = abs_ci_results[, 1], upper_abs = abs_ci_results[, 2])
        
        # Generazione del dataframe con i confronti a coppie sui valori assoluti
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
        
        # Confronti sui valori assoluti
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
        
        # Mostra gli istogrammi per ciascun trait
        if (show_histograms) {
                histograms <- lapply(1:nrow(effect_sizes), function(i) {
                        ggplot(data.frame(effect_dist = effect_sizes[i,]), aes(x = effect_dist)) +
                                geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
                                geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
                                labs(title = paste("Effect Size Distribution for", vars[i]), x = "Effect Size", y = "Frequency") +
                                theme_minimal()
                })
                
                # Stampa gli istogrammi
                for (hist in histograms) {
                        print(hist)
                }
        }
        
        # Mostra i confronti a coppie in base all'argomento
        if (show_pairwise_comparisons == TRUE) {
                # Visualizzazione dei confronti a coppie sui valori originali
                pairwise_plot <- ggplot(pairwise_comparisons, aes(x = Trait1, y = Trait2, fill = Significant_Difference)) +
                        geom_tile(color = "white") +
                        scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "grey")) +
                        labs(title = "Pairwise Significant Differences Between Traits (Original ES)", x = "Trait 1", y = "Trait 2") +
                        theme_minimal() +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                print(pairwise_plot)
                
                # Visualizzazione dei confronti a coppie sui valori assoluti
                abs_pairwise_plot <- ggplot(abs_pairwise_comparisons, aes(x = Trait1, y = Trait2, fill = Abs_Significant_Difference)) +
                        geom_tile(color = "white") +
                        scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "grey")) +
                        labs(title = "Pairwise Significant Differences Between Traits (Absolute ES)", x = "Trait 1", y = "Trait 2") +
                        theme_minimal() +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                print(abs_pairwise_plot)
                
        } else if (show_pairwise_comparisons == "PAIR") {
                # Visualizzazione dei confronti a coppie sui valori originali
                pairwise_plot <- ggplot(pairwise_comparisons, aes(x = Trait1, y = Trait2, fill = Significant_Difference)) +
                        geom_tile(color = "white") +
                        scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "grey")) +
                        labs(title = "Pairwise Significant Differences Between Traits (Original ES)", x = "Trait 1", y = "Trait 2") +
                        theme_minimal() +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
                print(pairwise_plot)
                
        } else if (show_pairwise_comparisons == "ABS") {
                # Visualizzazione dei confronti a coppie sui valori assoluti
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


results<-BootstrapEffectSizesOPTG(test_data,n_iter = 100,method="Cohens", GR= "group",show_histograms = F, show_pairwise_comparisons = F)
results
