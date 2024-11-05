library(dplyr)
library(ggplot2)
library(tidyr)
library(boot)


# Funzione per generare il dataset con NA casuali e nomi buffi presi dalla cultura pop
#si possono settare:
#n_samples = 100 (numero di campioni in ogni gruppo della popolazione)
#seed = (numero casuale)
#na_percentage = 0 (percentuale (da 0 a 1) di NA presenti del dataset)
#
#

generate_funny_dataset <- function(n_samples =60, seed = sample(1000000, 1), na_percentage = 0.1) {
        set.seed(seed)
        
        # Genera la colonna 'group' e variabili indipendenti con nomi buffi
        group <- rep(c("high", "low"), each = n_samples)
        mario <- c(rnorm(n_samples, mean = 10, sd = 6), rnorm(n_samples, mean = 8, sd = 2)) # Normal - Super Mario
        darwin <- c(rnorm(n_samples, mean = 8, sd = 10), rnorm(n_samples, mean = 8, sd = 11)) # Normal - Darwin but same mean
        elvis <- c(runif(n_samples, min = 5, max = 20), runif(n_samples, min = 4, max = 12))  # Uniform - Elvis Presley
        yoda <- c(rpois(n_samples, lambda = 10), rpois(n_samples, lambda = 5))                 # Poisson - Yoda (Star Wars)
        pacman <- c(rbinom(n_samples, size = 20, prob = 0.3), rbinom(n_samples, size = 20, prob = 0.5)) # Binomial - Pac-Man
        meme <- c(rexp(n_samples, rate = 0.5), rexp(n_samples, rate = 0.3))                   # Exponential - Meme culture
        
        # Combina in un data frame
        funny_dataset <- data.frame(group, mario, darwin, elvis, yoda, pacman, meme)
        
        # Introduce NA in modo casuale
        for (col in names(funny_dataset)[-1]) {  # Escludendo 'group'
                funny_dataset[[col]][sample(1:(2 * n_samples), size = round(na_percentage * 2 * n_samples))] <- NA
        }
        
        return(funny_dataset)
}


test_data <- generate_funny_dataset(n_samples = 50, na_percentage = 0.01, seed=123)
head(test_data)

###########################################################################################
# Test della funzione e significato dei parametri. Quelli con * sono obbligatori,
#gli altri non sono obbligatori e se non specificati utilizzano i valori indicati:
# * x= il dataset da usare
#Replacement = TRUE better to keep it as it is
#GR = "group" the name of the grouping variable (e.g. line, treatment, ecc)
#IN = 2 the first column where there are traits
#FIN = ncol(x) the last column where there are traits
#G1 = "high" the "tag" of individuals of one group (e.g. "ad libitum" or "diet"). 
#G2 = "low" the "tag" of individuals of the other group (e.g. "ad libitum" or "diet").
#n_iter = 1000
#conf_method = "percentile". can be also "Bca" "SD", or "bootstrap-t". these two however need to be tested and do not work perfectly
#SEED=123 change seed to get a different bootstrap
#method == "Cohen" or "Hedges"




################################################################
####Original Function
###############################################################

BootstrapEffectSizesG <- function(x, Replacement = TRUE, GR = "group", IN = 2, FIN = ncol(x), G1 = "high", G2 = "low", 
                                  n_iter = 1000, conf_method = "percentile", SEED = 123, method = "Hedges") {
        # Imposta il seed a livello della funzione principale
        set.seed(SEED)
        # Internal function to calculate effect sizes with NA handling for each variable
        NewPop <- function(x, Replacement, GR, IN, FIN, G1, G2) {
                
                pop1 <- x[x[[GR]] == G1, ]   # Filter for G1
                pop2 <- x[x[[GR]] == G2, ]   # Filter for G2
                
                effctsize <- numeric(FIN - IN + 1)  # Initialize a vector for effect sizes
                n1 <- numeric(FIN - IN + 1)
                n2 <- numeric(FIN - IN + 1)
                
                # Resampling for each variable independently
                for (i in IN:FIN) {
                        var_name <- names(x)[i]
                        
                        # Population 1 for the current variable
                        pop1_var <- pop1[[var_name]]
                        n1_var <- sum(!is.na(pop1_var))  # Effective size without NAs
                        resp1 <- sample(pop1_var[!is.na(pop1_var)], size = n1_var, replace = Replacement)
                        
                        # Population 2 for the current variable
                        pop2_var <- pop2[[var_name]]
                        n2_var <- sum(!is.na(pop2_var))  # Effective size without NAs
                        resp2 <- sample(pop2_var[!is.na(pop2_var)], size = n2_var, replace = Replacement)
                        
                        # Calculate means and standard deviations for each resampled variable
                        p1mean <- mean(resp1, na.rm = TRUE)
                        p2mean <- mean(resp2, na.rm = TRUE)
                        p1sd <- sd(resp1, na.rm = TRUE)
                        p2sd <- sd(resp2, na.rm = TRUE)
                        
                        # Calculate the effect size for the current variable
                        if (method == "Cohen") {
                                effctsize[i - IN + 1] <- (p1mean - p2mean) / sqrt(((n1_var - 1) * p1sd^2 + (n2_var - 1) * p2sd^2) / 
                                                                                          (n1_var + n2_var - 2))
                        } else if (method == "Hedges") {
                                effctsize[i - IN + 1] <- ((p1mean - p2mean) / sqrt(((n1_var - 1) * p1sd^2 + (n2_var - 1) * p2sd^2) / 
                                                                                           (n1_var + n2_var - 2))) * (1 - (3 / (4 * (n1_var + n2_var) - 9)))
                        }
                        # Save the effective sizes without NAs for each variable
                        n1[i - IN + 1] <- n1_var
                        n2[i - IN + 1] <- n2_var
                }
                
                return(list(effctsize = effctsize, n1 = n1, n2 = n2))
        }
        
        # Bootstrap distribution of effect sizes and sample sizes
        effect_size_boot <- replicate(n_iter, {
                result <- NewPop(x, Replacement, GR, IN, FIN, G1, G2)
                list(effctsize = result$effctsize, n1 = result$n1, n2 = result$n2)
        }, simplify = FALSE)
        
        # Extract effect size statistics
        effect_summary <- sapply(effect_size_boot, function(res) res$effctsize)
        mean_effect_summary <- rowMeans(effect_summary, na.rm = TRUE)
        
        ci_summary <- data.frame()
        histograms <- list()
        
        # Calculate confidence intervals
        for (i in 1:(FIN - IN + 1)) {
                effect_dist <- sapply(effect_size_boot, function(res) res$effctsize[i])
                
                ci <- c(NA, NA)
                
                if (conf_method == "SD") {
                        ci <- c(mean(effect_dist)-sd(effect_dist),mean(effect_dist)+sd(effect_dist)) #ci <- quantile(effect_dist, c(0.025, 0.975), na.rm = TRUE)
                } else if (conf_method == "percentile") {
                        ci <- quantile(effect_dist, c(0.025, 0.975), na.rm = TRUE)
                } else if (conf_method == "BCa") {
                        boot_obj <- boot(effect_dist, statistic = function(data, idx) mean(data[idx]), R = n_iter)
                        ci <- tryCatch({
                                boot_ci <- boot.ci(boot_obj, type = "bca")
                                if (!is.null(boot_ci$bca)) boot_ci$bca[4:5] else c(NA, NA)
                        }, error = function(e) c(NA, NA))
                } else if (conf_method == "bootstrap-t") {
                        boot_obj <- boot(effect_dist, statistic = function(data, idx) {
                                mean_val <- mean(data[idx])
                                se_val <- sd(data[idx]) / sqrt(length(data[idx]))
                                return(c(mean_val, se_val))
                        }, R = n_iter)
                        
                        ci <- tryCatch({
                                boot_ci <- boot.ci(boot_obj, type = "stud")
                                if (!is.null(boot_ci$student)) boot_ci$student[4:5] else c(NA, NA)
                        }, error = function(e) c(NA, NA))
                }
                
                # Mean sample sizes for the current variable
                mean_resampled_n1 <- mean(sapply(effect_size_boot, function(res) res$n1[i]), na.rm = TRUE)
                mean_resampled_n2 <- mean(sapply(effect_size_boot, function(res) res$n2[i]), na.rm = TRUE)
                
                # Summary for each variable
                var_name <- names(x)[IN + i - 1]
                ci_summary <- rbind(ci_summary, data.frame(
                        Variable = var_name,
                        mean = mean_effect_summary[i],
                        lower = ci[1],
                        upper = ci[2],
                        significant = ifelse(ci[1] > 0 | ci[2] < 0, TRUE, FALSE),
                        observed_n1 = sum(!is.na(x[x[[GR]] == G1, var_name])),
                        observed_n2 = sum(!is.na(x[x[[GR]] == G2, var_name])),
                        mean_resampled_n1 = mean_resampled_n1,
                        mean_resampled_n2 = mean_resampled_n2
                ), make.row.names = FALSE)
                
                # Create histogram of effect size distribution for the current variable
                hist_plot <- ggplot(data.frame(effect_dist), aes(x = effect_dist)) +
                        geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
                        geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
                        labs(title = paste("Effect Size Distribution for", var_name), x = "Effect Size", y = "Frequency") +
                        theme_minimal()
                
                histograms[[var_name]] <- hist_plot
        }
        
        # Plot of means and confidence intervals
        plot <- ggplot(ci_summary, aes(x = Variable, y = mean, ymin = lower, ymax = upper)) +
                geom_pointrange(color = "blue", size = 1) +
                geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
                labs(title = "Effect Sizes with Confidence Intervals", x = "Variables", y = "Effect Size") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        print(plot)
        # Print the method used for effect sizes, confidence method, and number of iterations
        cat("\nMethod used for Effect Sizes:", method, "\n")
        cat("Confidence Method:", conf_method, "\n")
        cat("Number of Iterations:", n_iter, "\n\n")
        # Return summary table and histograms
        return(list(summary = ci_summary, histograms = histograms))
}



##################################
##Versione alternativa semplificata con i loop 
####################################

BootstrapEffectSizesNEW <- function(x, Replacement = TRUE, GR = "group", IN = 2, FIN = ncol(x), G1 = "high", G2 = "low", 
                                    n_iter = 1000, conf_method = "percentile", SEED = 123, method = "Hedges") {
        set.seed(SEED)
        
        calc_effect_size <- function(pop1, pop2, method) {
                n1 <- length(pop1)
                n2 <- length(pop2)
                p1mean <- mean(pop1, na.rm = TRUE)
                p2mean <- mean(pop2, na.rm = TRUE)
                p1sd <- sd(pop1, na.rm = TRUE)
                p2sd <- sd(pop2, na.rm = TRUE)
                
                pooled_sd <- sqrt(((n1 - 1) * p1sd^2 + (n2 - 1) * p2sd^2) / (n1 + n2 - 2))
                effect <- (p1mean - p2mean) / pooled_sd
                
                if (method == "Hedges") {
                        effect <- effect * (1 - (3 / (4 * (n1 + n2) - 9)))
                }
                
                return(c(effect = effect, n1 = n1, n2 = n2))
        }
        
        vars <- names(x)[IN:FIN]
        pop1 <- x[x[[GR]] == G1, vars, drop = FALSE]
        pop2 <- x[x[[GR]] == G2, vars, drop = FALSE]
        
        boot_results <- replicate(n_iter, {
                sapply(vars, function(var) {
                        sample1 <- sample(pop1[[var]], sum(!is.na(pop1[[var]])), replace = Replacement)
                        sample2 <- sample(pop2[[var]], sum(!is.na(pop2[[var]])), replace = Replacement)
                        calc_effect_size(sample1, sample2, method)
                })
        }, simplify = "array")
        
        mean_effects <- rowMeans(boot_results["effect",, ], na.rm = TRUE)
        
        ci_summary <- lapply(vars, function(var) {
                x <- boot_results["effect", var, ]
                ci <- switch(conf_method,
                             "SD" = c(mean(x) - sd(x), mean(x) + sd(x)),
                             "percentile" = quantile(x, c(0.025, 0.975), na.rm = TRUE),
                             "BCa" = {
                                     boot_obj <- boot(x, function(d, i) mean(d[i]), R = n_iter)
                                     tryCatch(boot.ci(boot_obj, type = "bca")$bca[4:5], error = function(e) c(NA, NA))
                             },
                             "bootstrap-t" = {
                                     boot_obj <- boot(x, function(d, i) c(mean(d[i]), sd(d[i]) / sqrt(length(d[i]))), R = n_iter)
                                     tryCatch(boot.ci(boot_obj, type = "stud")$student[4:5], error = function(e) c(NA, NA))
                             })
                data.frame(
                        Variable = var,
                        mean = mean(x),
                        lower = ci[1],
                        upper = ci[2],
                        significant = ci[1] > 0 | ci[2] < 0,
                        observed_n1 = sum(!is.na(pop1[[var]])),
                        observed_n2 = sum(!is.na(pop2[[var]])),
                        mean_resampled_n1 = mean(boot_results["n1", var, ]),
                        mean_resampled_n2 = mean(boot_results["n2", var, ])
                )
        })
        
        ci_summary <- do.call(rbind, ci_summary)
        
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
        
        histograms <- lapply(vars, function(var) {
                ggplot(data.frame(effect_dist = boot_results["effect", var, ]), aes(x = effect_dist)) +
                        geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
                        geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
                        labs(title = paste("Effect Size Distribution for", var), x = "Effect Size", y = "Frequency") +
                        theme_minimal()
        })
        
        return(list(summary = ci_summary, histograms = histograms))
}





#################################################################################################################
#Versione ottimizzata per metterci poco tempo. Non è poi così diversa dall'originale
##################################################################################

BootstrapEffectSizesOPT <- function(x, Replacement = TRUE, GR = "group", IN = 2, FIN = ncol(x), G1 = "high", G2 = "low", 
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
        
        return(list(summary = ci_summary, histograms = histograms))
}




###########################################################################################
# Test della funzione e significato dei parametri. Quelli con * sono obbligatori,
#gli altri non sono obbligatori e se non specificati utilizzano i valori indicati:
# * x= il dataset da usare
#Replacement = TRUE better to keep it as it is
#GR = "group" the name of the grouping variable (e.g. line, treatment, ecc)
#IN = 2 the first column where there are traits
#FIN = ncol(x) the last column where there are traits
#G1 = "high" the "tag" of individuals of one group (e.g. "ad libitum" or "diet"). 
#G2 = "low" the "tag" of individuals of the other group (e.g. "ad libitum" or "diet").
#n_iter = 1000
#conf_method = "percentile". can be also "Bca" "SD", or "bootstrap-t". these two however need to be tested and do not work perfectly
#SEED=123 change seed to get a different bootstrap
#method == "Cohen" or "Hedges"


#######Updated version to compare effect sizes differences
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


results<-BootstrapEffectSizesOPT(test_data,n_iter = 100,method="Hedges", GR= "group")