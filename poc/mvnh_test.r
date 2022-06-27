library(MVNH)
library(tidyverse)


# Quick code (mostly written by Muyang) to test whether MVNH dissimilarity is 
# robust to violations of normality.  Answer: it is not...

sim.results = lapply(1:100, function(x){
  x1 = matrix(rlnorm(200,meanlog = 10, sdlog =2),ncol=2)
  x2 = matrix(rlnorm(200,meanlog = 12, sdlog =1),ncol=2)
  # simulate two variables from idd but only use the first variable for demonstration
  true.BD = MVNH_dissimilarity(x1,x2)
  true.BD = c(true.BD[[1]][2],true.BD[[2]][2],true.BD[[3]][2])
  norm.BD = MVNH_dissimilarity(log(x1),log(x2))
  norm.BD = c(norm.BD[[1]][2],norm.BD[[2]][2],norm.BD[[3]][2])
  # sqroot.BD = MVNH_dissimilarity(x1^0.5,x2^0.5)
  # sqroot.BD = c(sqroot.BD[[1]][2],sqroot.BD[[2]][2],sqroot.BD[[3]][2])
  # cube.BD = MVNH_dissimilarity(x1^3,x2^3)
  # cube.BD = c(cube.BD[[1]][2],cube.BD[[2]][2],cube.BD[[3]][2])
  return(c(true.BD,norm.BD))
})
sim.results = as.data.frame(t(matrix(unlist(sim.results),nrow = 6)))
names(sim.results) = paste0(rep(c("original (lognormal)","log(orginal)"),each=3),
                            c("_BD","_MD","_DR"))
par(mfrow=c(1,3))
boxplot(sim.results[,grep("_BD",names(sim.results))],ylab="BD",xlab="Transformation",main="n=100")
boxplot(sim.results[,grep("_MD",names(sim.results))],ylab="MD",xlab="Transformation",main="n=100")
boxplot(sim.results[,grep("_DR",names(sim.results))],ylab="DR",xlab="Transformation",main="n=100")

mean(sim.results$`original (lognormal)_MD`)/mean(sim.results$`original (lognormal)_DR`)
mean(sim.results$`log(orginal)_MD`)/mean(sim.results$`log(orginal)_DR`)
