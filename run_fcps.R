#!/usr/bin/env/R

## Omnibenchmark-izes Marek Gagolewski's https://github.com/gagolews/clustering-results-v1/blob/eae7cc00e1f62f93bd1c3dc2ce112fda61e57b58/.devel/do_benchmark_fcps_aux.R

## Takes the true number of clusters into account and outputs a 2D matrix with as many columns as ks tested,
## being true number of clusters `k` and tested range `k plusminus 2`


library(argparse)
library(FCPS)
## library(R.utils)

parser <- ArgumentParser(description="FCPS caller")


parser$add_argument('--data.matrix',
                    type="character",
                    help='gz-compressed textfile containing the comma-separated data to be clustered.')
parser$add_argument('--data.true_labels',
                    type="character",
                    help='gz-compressed textfile with the true labels; used to select a range of ks.')
parser$add_argument('--seed',
                    type="integer",
                    help='Random seed',
                    default = 819797,
                    dest = 'seed')
parser$add_argument("--output_dir", "-o", dest="output_dir", type="character",
                    help="output directory where files will be saved", default=getwd())
parser$add_argument("--name", "-n", dest="name", type="character", help="name of the dataset")
parser$add_argument("--method", "-m", dest="method", type="character", help="method")

args <- parser$parse_args()


# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))


VALID_METHODS <- list(
    # Affinity propagation (Apclustering) - does not allow k
    # DBSCAN - does not allow k
    # DensityPeakClustering - does not allow k
    # MarkovClustering- does not allow k
    # MSTclustering - does not allow k
    # OPTICSclustering - does not allow k
    # PenalizedRegressionBasedClustering - does not allow k
    # pdfClustering - does not allow k
    # QTclustering - does not allow k
    # SharedNearestNeighborClustering - does not allow k
    # SubspaceClustering - does not allow k

    # AutomaticProjectionBasedClustering - finds a nonlinear projection (different class)
    #list(HierarchicalClustering, "Sparse") # does auto feature selection -> skip (separate class)
    # ProjectionPursuitClustering
    # RobustTrimmedClustering
    # SOMclustering
    # Spectrum - spectral clustering
    # SubspaceClustering - subspaces
    # TandemClustering - combines k-means and PCA

    # CrossEntropyClustering - FCPS documentation mentions that the algorithm is not stable
    # DatabionicSwarmClustering - current implementation is not efficient for N>4000
    # Agglomerative Nesting (AGNES) - this is the ordinary Hierarchical Clustering

    # ModelBasedClustering, MoGclustering - mixture of Gaussians - see scikit-learn

    FCPS_AdaptiveDensityPeak=list(ADPclustering), # default params
    FCPS_Minimax=list(HierarchicalClustering, "Minimax"),  # protoclust::protoclust Minimax Linkage; no params
    FCPS_MinEnergy=list(HierarchicalClustering, "MinEnergy"),  # energy::energy.hclust Hierarchical Clustering by Minimum (Energy) E-distance; MinimalEnergyClustering no params
    FCPS_HDBSCAN_4=list(HierarchicalClustering, "HDBSCAN"), # dbscan::hdbscan HierarchicalDBSCAN minPts=4
    FCPS_HDBSCAN_2=list(HierarchicalClustering, "HDBSCAN", minPts=2),
    FCPS_HDBSCAN_8=list(HierarchicalClustering, "HDBSCAN", minPts=8),
    FCPS_Diana=list(DivisiveAnalysisClustering), # cluster::diana DIvisive ANAlysis Clustering
    FCPS_Fanny=list(FannyClustering, maxit=2000), # cluster::fanny Fuzzy Analysis Clustering
    FCPS_Hardcl=list(HCLclustering), # cclust::cclust(method="hardcl") On-line Update (Hard Competitive learning convex clustering) method 
    FCPS_Softcl=list(NeuralGasClustering), # cclust::cclust(method="neuralgas")  Neural Gas (Soft Competitive learning) 
    
    FCPS_Clara=list(LargeApplicationClustering, Standardization=FALSE,Random=TRUE), # cluster::clara Clustering Large Applications - based on Partitioning Around Medoids on subsets; 
    FCPS_PAM=list(PAMclustering) #  cluster::pam Partitioning Around Medoids (PAM)
)


load_labels <- function(data_file) {
    (fd <- read.table(gzfile(data_file), header = FALSE)$V1)
}

load_dataset <- function(data_file) {
    (fd <- read.table(gzfile(data_file), header = FALSE))
}


pin_seed <- function(fun, args, seed) {
    set.seed(seed)
    print(fun)
    eval(as.call(c(fun, args)))
    cat("after eval() in pin_seed()\n")
}

do_fcps <- function(data, Ks, method, seed) {
    if (!method %in% names(VALID_METHODS))
        stop('Not a valid method')
    
    d <- as.matrix(dist(data))
    data <- as.matrix(data)
    
    res <- list()
    case <- VALID_METHODS[[method]]
    fun <- case[[1]]

    last_k <- 0
    for (i in 1:length(Ks)) {
        cat('current Ks index is ', i, 'and k is ', k, '\n')
        k <- Ks[i]
        
        if ("DataOrDistances" %in% names(formals(fun)))
            args <- list(DataOrDistances=d)
        else
            args <- list(Data=data)        
        
        args <- c(args, ClusterNo=k, case[-1])
        print(lapply(args, head))

        y_pred <- as.integer(pin_seed(fun, args, seed)[['Cls']])

        print(table(y_pred))
        if (min(y_pred) > 0 && max(y_pred) == k) {
            res[[paste('k ', k, 'i ',  i)]] <- y_pred
        }
        else {
            ## error means all are assigned to the same cluster
            res[[paste('k ', k, 'i ', i)]] <- rep(k, length(y_pred))
        }

        ## double check repeated ks are consistent
        if (k == last_k) {
            stopifnot(identical(res[[i]], res[[(i-1)]]))
            cat('   repeated k results are consistent\n')
            
        }
        last_k <- k
       
    }
    
    return(do.call('cbind.data.frame', res))
}

truth <- load_labels(args[['data.true_labels']])

k <- max(truth) # true number of clusters
Ks <- c(k-2, k-1, k, k+1, k+2) # ks tested, including the true number

cat('original Ks grid: ', Ks, '\n')

Ks[Ks < 2] <- 2 ## but we never run k < 2; those are replaced by (extra) k=2 runs (not to skip the calculation)

cat('refined Ks grid:  ', Ks, '\n')


dataset <- load_dataset(args[['data.matrix']])
head(dataset)
dim(dataset)

res <- do_fcps(data = dataset, 
               method = args[['method']], 
               Ks = Ks, seed = args$seed)

colnames(res) <- paste0('k=', Ks)
    
gz <- gzfile(file.path(args[['output_dir']], paste0(args[['name']], "_ks_range.labels.gz")), "w")
write.table(file = gz, res, col.names = TRUE, row.names = FALSE, sep = ",")
close(gz)
