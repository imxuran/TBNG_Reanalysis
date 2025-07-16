### Conduct cPCA with different alphas


alphas <- 10^seq(0, 1, length.out = 7) 
#alphas <- 1:10

lapply(alphas, cPCA_plot)