# Background

Here we go over a few basic concepts in Population Genetics and in Statistics that are useful for grasping the subleties of CTGA.

Some popular R packages for population genetics
Population genetics in R
http://grunwaldlab.github.io/Population_Genetics_in_R/

PopGenome: An Efficient Swiss Army Knife for Population Genomic Analyses
https://cran.r-project.org/web/packages/PopGenome/

### Genetic Drift
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

### Metrics of nucleotide diversity

There are numerous metrics for DNA variability. The two main ones are Watterson’s estimate ([theta](https://en.wikipedia.org/wiki/Watterson_estimator)) and Tajima’s diversity ([pi](https://en.wikipedia.org/wiki/Nucleotide_diversity)). Watterson's metrics depends on the number of SNPs per sequence length. Tajima’s diversity estimate is the average number of allele differences between all pairs of individuals. In all cases, the estimator must be divided by the number of base pairs analyzed.

Under a neutral model (i.e., only drift, panmixia and constant effective size), both estimators should be the same.
  
#### TOY EXERCISE
Here is a list of 7 aligned sequences

    AAATTTTCCGGCA
    .............
    ..T.........G
    .............
    ......A......
    ......A......
    ..T..........
    
-	How 	many SNPs?
-	Compute theta and pi diversity

#### Real example

The attached file [MC1R_pigs.fasta](https://github.com/miguelperezenciso/CTGAcourse/blob/master/Chapter1_Background/MC1R_Pigs_aligned.fasta) contains a list of pig MC1R gene sequences (Fang et al. 2018, Contrasting Mode of Evolution at a Coat Color Locus in Wild and Domestic Pigs, https://doi.org/10.1371/journal.pgen.1000341). A very popular software to analyze sequence data in a small scale is DNASp (http://www.ub.edu/dnasp/). It runs only in windows.

- Install DNAsp the program and upload fasta file.
- Compute diversity estimates for the whole set of samples and separately by continent (Europe vs. Asia). 
- Compare both Tajima's and Watterson's estimators.

# Breeding

- Write an R function that returns mean of populations given allele effects and frequencies
- Same but returning mean of crosses between two lines

# Simulation
Simulation is a key tool in many quantitative genetics applications. Give a thought at how would you design a program to simulate a F2 cross between two inbred lines.
