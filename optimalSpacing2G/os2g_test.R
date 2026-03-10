library(secrdesign)
source('optimalSpacing2G/optimalSpacing2G.R')
grid <- make.grid(7,7,50)
os1 = optimalSpacing(D = 5, traps = grid, detectpar = list(lambda0 = 0.2, sigma = 40), 
               noccasions = 5, plt = FALSE)
os1$rotRSE$optimum.spacing

os2 = optimalSpacing2G(D1 = 5, D2 = 5,
                 traps0 = grid,
                 detectpar1 = list(lambda0 = 0.2, sigma = 40),
                 detectpar2 = list(lambda0 = 0.2, sigma = 40),
                 noccasions = 5,
                 criterion = c("sum_min"), # all_min
                 spacing_m = seq(0,120,1))

os2$optimum.spacing
xx=os2$values
xx[1,]
