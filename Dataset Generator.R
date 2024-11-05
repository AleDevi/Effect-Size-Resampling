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