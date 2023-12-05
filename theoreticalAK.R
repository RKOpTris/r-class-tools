# predict accuracy/kappa when classes are balanced

theoreticalAK <- function(accuracy = NULL, kappa = NULL, n_class){
  if(!is.null(kappa)){
    accuracy <- ((1 - (1 / n_class)) * kappa) + (1 / n_class)
    return(accuracy)
  } else if(!is.null(accuracy)){
    kappa <- (accuracy - (1 / n_class)) / (1 - 1 / n_class)
    return(kappa)
  }
}

# theoreticalAK(kappa = 0, n_class = 2)
# theoreticalAK(accuracy = c(0.8, 0.9), n_class = 2)
