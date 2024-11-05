# Carica i pacchetti necessari
library(ggplot2)
library(dplyr)

# Definisci le funzioni originali se non già definite
# BootstrapEffectSizesG, BootstrapEffectSizesNEW, BootstrapEffectSizesOPT

# Parametri per il benchmark
sample_sizes <- seq(from = 0, to = 250, by = 100)  # Campioni di grandezza crescente
num_traits <- c(5)  # Numero di tratti da considerare
iterations <- c(500, 1000)  # Numero di iterazioni

# Crea un data frame per memorizzare i risultati
benchmark_results <- data.frame()

# Lista delle funzioni con i loro nomi
functions <- list(
        "BootstrapEffectSizesG" = BootstrapEffectSizesG,
        "BootstrapEffectSizesNEW" = BootstrapEffectSizesNEW,
        "BootstrapEffectSizesOPT" = BootstrapEffectSizesOPT
)

# Ciclo per le varie dimensioni del campione, numero di tratti e iterazioni
for (size in sample_sizes) {
        for (traits in num_traits) {
                for (iter in iterations) {
                        # Crea un dataset di test con dimensioni variabili
                        test_data <- data.frame(
                                group = rep(c("high", "low"), each = size / 2),
                                matrix(rnorm(size * traits), ncol = traits)
                        )
                        colnames(test_data)[-1] <- paste0("Trait_", 1:traits)
                        
                        # Valuta ciascuna funzione e registra i tempi
                        for (func_name in names(functions)) {
                                # Seleziona la funzione specifica
                                fn <- functions[[func_name]]
                                
                                # Calcola il tempo di esecuzione
                                exec_time <- system.time(
                                        fn(test_data, GR = "group", IN = 2, FIN = 1 + traits, n_iter = iter)
                                )[3]
                                
                                # Salva i risultati
                                benchmark_results <- rbind(benchmark_results, data.frame(
                                        Function = func_name,
                                        Sample_Size = size,
                                        Traits = traits,
                                        Iterations = iter,
                                        Time = exec_time
                                ))
                        }
                }
        }
}
# Visualizzazione dei risultati
ggplot(benchmark_results, aes(x = Sample_Size, y = Time, color = Function)) +
        geom_line() +
        geom_point() +
        facet_grid(Traits ~ Iterations, scales = "free_y") +
        labs(
                title = "Benchmark delle Funzioni BootstrapEffectSizes",
                x = "Numero di Campioni",
                y = "Tempo di Esecuzione (secondi)"
        ) +
        theme_minimal() +
        theme(legend.position = "top")

# Visualizzazione dei risultati in un unico grafico
ggplot(benchmark_results, aes(x = Iterations, y = Time, color = Function, shape = as.factor(Traits), linetype = as.factor(Sample_Size))) +
        geom_line() +
        geom_point(size = 2) +
        labs(
                title = "Benchmark delle Funzioni BootstrapEffectSizes",
                x = "Numero di Campioni",
                y = "Tempo di Esecuzione (secondi)",
                color = "Funzione",
                shape = "Numero di Tratti",
                linetype = "Iterazioni"
        ) +
        theme_minimal() +
        theme(legend.position = "top")
