result.acc <- rep(0,25)
result.ess <- rep(0,25)
setwd("Q:/high-dimension-sircovid-MRes/draft/vaccine_delay_fits_gradient")
library(coda)


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

df <- data.frame(x = rep(1:25,1)/100)

################ L 1 ################
L <- 1
for (i in 1:25){
  x <- read.csv(sprintf('L_%s_1/outputs/samples_L_%s_epsilon_%s.csv',L,L,i/100))
  x <- as.matrix(x[250:1000,2:29])
  for (j in 2:10){
    temp <- read.csv(sprintf('L_%s_%s/outputs/samples_L_%s_epsilon_%s.csv',L,j,L,i/100))
    temp <- as.matrix(temp[250:1000,2:29])
    x <- rbind(x,temp)
  }
  result.acc[i] <- acc_rate(x)
  result.ess[i] <- sum(effectiveSize(x))
}
df$acc.1 <- result.acc
df$ess.1 <- result.ess


################ L 2 ################
L <- 2
for (i in 1:25){
  x <- read.csv(sprintf('L_%s_1/outputs/samples_L_%s_epsilon_%s.csv',L,L,i/100))
  x <- as.matrix(x[250:1000,2:29])
  for (j in 2:10){
    temp <- read.csv(sprintf('L_%s_%s/outputs/samples_L_%s_epsilon_%s.csv',L,j,L,i/100))
    temp <- as.matrix(temp[250:1000,2:29])
    x <- rbind(x,temp)
  }
  result.acc[i] <- acc_rate(x)
  result.ess[i] <- sum(effectiveSize(x))
}
df$acc.2 <- result.acc
df$ess.2 <- result.ess

################ L 3 ################
L <- 3
for (i in 1:25){
  x <- read.csv(sprintf('L_%s_1/outputs/samples_L_%s_epsilon_%s.csv',L,L,i/100))
  x <- as.matrix(x[250:1000,2:29])
  for (j in 2:10){
    temp <- read.csv(sprintf('L_%s_%s/outputs/samples_L_%s_epsilon_%s.csv',L,j,L,i/100))
    temp <- as.matrix(temp[250:1000,2:29])
    x <- rbind(x,temp)
  }
  result.acc[i] <- acc_rate(x)
  result.ess[i] <- sum(effectiveSize(x))
}
df$acc.3 <- result.acc
df$ess.3 <- result.ess

p1 <- ggplot(df,aes(x = x)) + 
  geom_line(aes(y = acc.1,color = "L = 1"),lwd=1.2)+
  geom_line(aes(y = acc.2,color = "L = 2"),lwd=1.2)+
  geom_line(aes(y = acc.3,color = "L = 3"),lwd=1.2)+
  scale_fill_manual("", values = c("L = 1" = "blue", "L = 2" = "violet", "L = 3" = "red"))+
  ggtitle("Acceptance Rate of Hamiltonian Monte Carlo with L=1,2,3")+
  labs(x = "h" , y =  "acceptance rate",
       color = "Legend")
p1

p2 <- ggplot(df,aes(x = x)) + 
  geom_line(aes(y = ess.1,color = "L = 1"),lwd=1.2)+
  geom_line(aes(y = ess.2,color = "L = 2"),lwd=1.2)+
  geom_line(aes(y = ess.3,color = "L = 3"),lwd=1.2)+
  scale_fill_manual("", values = c("L = 1" = "blue", "L = 2" = "violet", "L = 3" = "red"))+
  ggtitle("Effective Sample Size of Hamiltonian Monte Carlo with L=1,2,3")+
  labs(x = "h" , y =  "effective sample size",
       color = "Legend")
p2


t1 <- 10.7/25
t2 <- (8.22+7.76)/25
t3 <-(12.4+8.24)/25
t = c(t1,t2,t3)

df$ess_time.1 <- df$ess.1/t1
df$ess_time.2 <- df$ess.2/t2
df$ess_time.3 <- df$ess.3/t3

p3 <- ggplot(df,aes(x = x)) + 
  geom_line(aes(y = ess_time.1/10,color = "L = 1"),lwd=1.2)+
  geom_line(aes(y = ess_time.2/10,color = "L = 2"),lwd=1.2)+
  geom_line(aes(y = ess_time.3/10,color = "L = 3"),lwd=1.2)+
  scale_fill_manual("", values = c("L = 1" = "blue", "L = 2" = "violet", "L = 3" = "red"))+
  ggtitle("Effective Sample Size per hour of Hamiltonian Monte Carlo with L=1,2,3")+
  labs(x = "h" , y =  "effective sample size per hour",
       color = "Legend")
p3

p <- ggplot(df2,aes(x = x, y =t, fill = x))+geom_bar(stat="identity")+ggtitle("running time per 1000 iterations")+labs(x = 'L',y='running time [hour]')
p
