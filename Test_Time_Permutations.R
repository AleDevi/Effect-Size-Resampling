library(dplyr)
library(ggplot2)
library(tidyr)
library(boot)

#Definizione delle tre funzioni (originale, alternativa e ottimizzata)

# Funzione originale
BootstrapEffectSizesG <- function(x, Replacement = TRUE, GR = "group", IN = 2, FIN = ncol(x), G1 = "high", G2 = "low", 
                                  n_iter = 1000, conf_method = "percentile", SEED = 123, method = "Hedges") {
        # ... [inserisci qui il codice completo della funzione originale] se non è gia in memoria...
}

# Funzione alternativa
BootstrapEffectSizesNEW <- function(x, Replacement = TRUE, GR = "group", IN = 2, FIN = ncol(x), G1 = "high", G2 = "low", 
                                    n_iter = 1000, conf_method = "percentile", SEED = 123, method = "Hedges") {
        # ... [inserisci qui il codice completo della funzione alternativa] se non è gia in memoria...
}

# Funzione ottimizzata
BootstrapEffectSizesOPT <- function(x, Replacement = TRUE, GR = "group", IN = 2, FIN = ncol(x), G1 = "high", G2 = "low", 
                                    n_iter = 1000, conf_method = "percentile", SEED = 123, method = "Hedges") {
        # ... [inserisci qui il codice completo della funzione ottimizzata] se non è gia in memoria...
}

# Funzione per testare le prestazioni
test_performance <- function(data, iterations) {
        results <- data.frame(
                expr = c("NEW", "ORIGINAL", "OPT"),
                mean = c(
                        system.time(BootstrapEffectSizesNEW(data, n_iter = iterations))[3],
                        system.time(BootstrapEffectSizesG(data, n_iter = iterations))[3],
                        system.time(BootstrapEffectSizesOPT(data, n_iter = iterations))[3]
                ),
                iterations = iterations
        )
        return(results)
}

# Creazione di un dataset di esempio
set.seed(123)
test_data <- data.frame(
        group = rep(c("high", "low"), each = 30),
        var1 = rnorm(60),
        var2 = rnorm(60),
        var3 = rnorm(60)
)

# Test con diverse numerosità di ricampionamento
iterations_to_test <- c(100, 500)
#iterations_to_test <- seq(from = 25000, to = 25000, by = 5000)
performance_results <- lapply(iterations_to_test, function(iter) {
        test_performance(test_data, iter)
})

all_results <- do.call(rbind, performance_results)

all_results <- all_results %>%
        mutate(
                seconds = mean,
                method = factor(expr, levels = c("NEW", "ORIGINAL", "OPT"), 
                                labels = c("Originale", "Alternativa", "Ottimizzata"))
        )
# Creazione del grafico

ggplot(all_results, aes(x = iterations, y = seconds, color = method)) +
        geom_line() +
        geom_point() +
        labs(title = "Confronto delle prestazioni delle funzioni Bootstrap",
             x = "Numero di iterazioni",
             y = "Tempo di esecuzione medio (secondi)",
             color = "Metodo") +
        theme_minimal() +
        scale_x_continuous(breaks = iterations_to_test) +
        scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

# Salvataggio del grafico

#ggsave("bootstrap_performance_comparison.png", width = 10, height = 6)
print(all_results %>%
              select(method, iterations, seconds) %>%
              arrange(iterations, method) %>%
              mutate(seconds = round(seconds, 3)) %>%
              pivot_wider(names_from = method, values_from = seconds) %>%
              arrange(iterations))
