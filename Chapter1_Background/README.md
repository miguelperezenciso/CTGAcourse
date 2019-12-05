# Background

Here we go over a few basic concepts in Population Genetics and in Statistics that are useful for grasping the subleties of CTGA.

Some popular R packages for population genetics
Population genetics in R
http://grunwaldlab.github.io/Population_Genetics_in_R/

PopGenome: An Efficient Swiss Army Knife for Population Genomic Analyses
https://cran.r-project.org/web/packages/PopGenome/

Simple genetic drift program

``
  #----- simulates drift
  # N is population size and p is initial allele frequency
  drift <- function(N,p) {
      g=as.numeric()
      t=0
      while(p>0 & p<1) {
      genotypes <- rbinom(N,1,p)
      p<-mean(genotypes)
      t=t+1
      # the last column of each row contains current allele frequency
      genotypes=cbind(t(genotypes),p)
      g=rbind(g,genotypes)
    }
    return(g[,N+1])
  }
``

EXERCISE: Guess what the program does and why
HINT: look for help with ‘help(command)’ or ‘? command’
