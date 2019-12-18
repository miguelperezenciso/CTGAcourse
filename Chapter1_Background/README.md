![](https://github.com/miguelperezenciso/CTGAcourse/blob/master/Chapter1_Background/Screenshot%20from%202019-12-18%2013-44-59.png)
# Background
In this chapter, we go over basic concepts in Population Genetics, Statistics and Bionformatics that are useful for grasping the subleties of CTGA. This section is covered by three ppt files:

- [Population Genetics and Statistics concepts](https://github.com/miguelperezenciso/CTGAcourse/blob/master/Chapter1_Background/CTGA_Chap1a_Background_pptx.pdf)
- [Breeding Basics](https://github.com/miguelperezenciso/CTGAcourse/blob/master/Chapter1_Background/CTGA_Chap1b_Breeding_pptx.pdf)
- [Some Bioinformatics Applications](https://github.com/miguelperezenciso/CTGAcourse/blob/master/Chapter1_Background/CTGA_Chap1c_Bioinformatics_pptx.pdf)

# 1. Population Genetics
Ultimately, the goal is to identify genetic and phenotypic variabilities. The sources of genetic variation are:

- Mutation: the ultimate source of all variability, due to replication errors. Note that only germline mutations are relevant for our purposes. Mutation is a rare process, from where it follows that SNP frequencies are also low.
- Recombination: Much has been debated on the evolutionary advantages of recombination; it allows shuffling extant allele combinations and improve adaptation to changing environments. NOTE: GWAS only apply in the presence of recombination. GWAS and Y or mitochondrial genomes make no sense.
- Genetic drift: It refers to the random sampling of haplotypes transmitted to offspring. Its strength is governed by the famous concept of **effective population size**.
- Selection: In this context, we are mainly interested in artificial selection, also called **directional selection**. There are many variants of selection: balancing, background, etc. Selection affects only a small subset of loci, whereas drift influences the whole genome.
- Migration refers to individuals moving from one population to another, and causes **admixture** or **genetic structure**. You can visualize structure with a principal component analyses (PCA).

Some popular R packages for population genetics:

http://grunwaldlab.github.io/Population_Genetics_in_R/

PopGenome: An Efficient Swiss Army Knife for Population Genomic Analyses
https://cran.r-project.org/web/packages/PopGenome/


### Metrics of nucleotide diversity
Broadly, mutation generates new variability whereas drift erodes it. A basic question is then measuring genetic variability. There are numerous metrics for DNA variability. The two main ones are Watterson’s estimate ([theta](https://en.wikipedia.org/wiki/Watterson_estimator)) and Tajima’s diversity ([pi](https://en.wikipedia.org/wiki/Nucleotide_diversity)). Watterson's metrics depends on the number of SNPs per sequence length. Tajima’s diversity estimate is the average number of allele differences between all pairs of individuals. In all cases, the estimator must be divided by the number of base pairs analyzed. Under a neutral model (i.e., only drift, panmixia and constant effective size), both estimators should be the same. Differences between both estimators appear because of selection or demographic events (bottlenecks, admixture, ...).
  
#### Toy exercise
Here is a list of 7 aligned sequences

    AAATTTTCCGGCA
    .............
    ..T.........G
    .............
    ......A......
    ......A......
    ..T..........
    
-	How many SNPs?
-	Compute theta and pi diversity

#### Real example
The attached file [MC1R_pigs.fasta](https://github.com/miguelperezenciso/CTGAcourse/blob/master/Chapter1_Background/MC1R_Pigs_aligned.fasta) contains a list of pig MC1R gene sequences (Fang et al. 2018, Contrasting Mode of Evolution at a Coat Color Locus in Wild and Domestic Pigs, https://doi.org/10.1371/journal.pgen.1000341). A very popular software to analyze sequence data in a small scale is DNASp (http://www.ub.edu/dnasp/). **It runs only in windows**.

- Install DNAsp the program and upload fasta file.
- Compute diversity estimates for the whole set of samples and separately by continent (Europe vs. Asia). 
- Compare both Tajima's and Watterson's estimators.

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


EXERCISE
- Guess what the program does and why.

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

### Admixture
The degree of differentiation between populations is usually measured with the Fst statistics. Sewall Wright defined Fst as the correlation between gametes chosen randomly from within the same subpopulation relative to the entire population. (https://www.nature.com/articles/nrg2611) 

EXERCISE:
- Write an R function to compute genotype frequencies of an admixed population. Parameters needed are allele frequencies in each population (p1, p2), and fraction of individuals from population 1 in admixed population (m1, logically m2=1-m1). We assume Hardy-Weinberg equilibrium within each founder population.

## 2. Statistics: "In God we trust, all others bring data"*
* Found in Found in Tibshirani et al., attributed to both Deming and Heyden

The term 'Statistics' and 'State' share the same etimology, as Statistics was initially (end 18th century) a means for the State to have censuses that make it possible to collect taxes as efficiently as possible. The real importance of Statistics in science lies in its role as quantifying uncertainty. In modern times, the concepts of **inference** and **model** are central to this goal. You should not believe that scientists are necessarily rational, well behaved people. Statistics is not a coherent and unified framework, many different schools (frequentist, Bayesian, non parametric) coexist and many angry disputes have been witnessed.

## Statistical Model
* A model is an abstraction of reality, it never exists a true model.
* 'The value of a model is that it often suggests a simple summary of data in terms of the major systematic effects together with a summary of the nature and amount of unexplained variation' (McCullagh and Nelder, 1983).
* Two desirable characteristics of a model are parsimony and goodness of fit. Here, the key is to find an optimum goodness of fit with the minimum number of parameters that maximizes the predictive ability of the model.
* More importantly a model in Quantitative Genetics is also a definition of a trait. For instance, the image shows a QTL profile with the same phenotype corrected by either age or weight ([Pérez-Enciso et al.](https://academic.oup.com/jas/article/78/10/2525/4670866)).

![](https://github.com/miguelperezenciso/CTGAcourse/blob/master/Chapter1_Background/Screenshot%20from%202019-12-18%2016-25-05.png)

## Multiple testing


## 3. Breeding

EXERCISES
- Write an R function that returns mean of populations given allele effects and frequencies
- Same but returning mean of crosses between two lines

# Simulation
Simulation is a key tool in many quantitative genetics applications. Give a thought at how would you design a program to simulate a F2 cross between two inbred lines.

## 4. Bioinformatics

