# Background

Here we go over a few basic concepts in Population Genetics and in Statistics that are useful for grasping the subleties of CTGA.

Some popular R packages for population genetics
Population genetics in R
http://grunwaldlab.github.io/Population_Genetics_in_R/

PopGenome: An Efficient Swiss Army Knife for Population Genomic Analyses
https://cran.r-project.org/web/packages/PopGenome/

Simple genetic drift program
 
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


EXERCISE: Guess what the program does and why

Running the drift R function

    # num of individuals
    N=20
    # initial allele frequency
    p0=.2
    f=drift(N,p0)
    # list of colors for plotting each replicate
    color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    plot(f, type='l', ylim=c(0,1), xlim=c(1,30), xlab='generation', ylab='f')
    for (rep in seq(100)) {
        f = drift(N, p0)
        lines(f, col=sample(color,1))
    }
    
EXERCISES 
- try different sizes (N) an dinitial frequencies (p)
- What is the probability of fixation of allele 1 ?

  
