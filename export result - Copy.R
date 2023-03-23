setwd("Q:/high-dimension-sircovid-MRes/draft/vaccine_delay_fits_gradient")
library(coda)
library(ggplot2)
acc_rate <- function(x){
  num <- dim(x)[1]-1
  count <- 0
  for (i in 1:num){
    if (x[i+1]==x[i]){
      count<-count+1
    }
  }
  return((num-count)/num)
}

x <- read.csv(sprintf('L_%s_1/outputs/samples_L_%s_epsilon_%s.csv',1,1,1/100))
x <- as.matrix(x[250:1000,2:29])
df <- data.frame(epsilon = 0.01, ess = sum(effectiveSize(x)), acc = acc_rate(x))

################ L 2 ################
L <- 3
for (i in 1:25){
  
  for (j in 1:10){
    
    x <- read.csv(sprintf('L_%s_%s/outputs/samples_L_%s_epsilon_%s.csv',L,j,L,i/100))
    x <- as.matrix(x[250:1000,2:29])
    df <- rbind(df, c(epsilon = i/100, ess = sum(effectiveSize(x)), acc = acc_rate(x)))
  }
}
df$epsilon <- as.factor(df$epsilon)
df <- df[2:251,]
p <-  ggplot(df, aes(x = epsilon, y=acc))+
  geom_boxplot(fill = "red")+
  ggtitle("Acceptance Rate of Hamiltonian Monte Carlo with L=3")+
  labs(x = "h" , y =  "acceptance rate")
p

