setwd("~/2024_hvCpG/")

## Load algorithm
source("05_hvCpGalgorithm/hvCpG_algorithm_detection_v5batches.R")

## Load hvCpG of Maria and matching mQTL
cistrans_GoDMC_hvCpG_matched_control <- read.table("03_prepDatasetsMaria/cistrans_GoDMC_hvCpG_matched_control.txt", header = T)

## Names of all CpG in the dataset:                                                                                                                                                   
cpg_names_all <<- rhdf5::h5read(paste0("/home/alice/arraysh5/all_matrix_noscale.h5"), "cpg_names")
length(cpg_names_all); head(cpg_names_all)
## [1] 394240
## [1] "cg00000029" "cg00000108" "cg00000109" "cg00000165" "cg00000236"
## [6] "cg00000289"

## In which positions are they?
cpg_names_all[cpg_names_all %in% cistrans_GoDMC_hvCpG_matched_control$hvCpG_name] %>% length
cpg_names_all[cpg_names_all %in% cistrans_GoDMC_hvCpG_matched_control$controlCpG_name]  %>% length

sub_cistrans_GoDMC_hvCpG_matched_control <- cistrans_GoDMC_hvCpG_matched_control[
  cistrans_GoDMC_hvCpG_matched_control$hvCpG_name %in% cpg_names_all &
    cistrans_GoDMC_hvCpG_matched_control$controlCpG_name %in% cpg_names_all,]

## Positions of targets in cpg_names:
pos2check <- c(match(sub_cistrans_GoDMC_hvCpG_matched_control$hvCpG_name, cpg_names_all),
               match(sub_cistrans_GoDMC_hvCpG_matched_control$controlCpG_name, cpg_names_all))

system.time(runAndSave(analysis = "arrayAlgov5", cpgPos_vec = head(pos2check,1000),
                       dataDir = "/home/alice/arraysh5",
                       resultDir="05_hvCpGalgorithm/resultsDir/MariaArrays/",
                       p0=0.80, p1=0.65, overwrite = TRUE,
                       NCORES=30, batch_size = 1000))

################
## Benchmarking:
library(bench)

# Wrap my function call
run_bench <- function(NCORES, batch_size) {
  runAndSave(
    analysis    = "arrayAlgov5",
    cpgPos_vec  = head(pos2check, 100000),
    dataDir     = "/home/alice/arraysh5",
    resultDir   = "05_hvCpGalgorithm/resultsDir/MariaArrays/",
    p0          = 0.80,
    p1          = 0.65,
    overwrite   = TRUE,
    skipsave    = TRUE,
    NCORES      = NCORES,
    batch_size  = batch_size
  )
}

# Use bench::press to test all combinations
results <- bench::press(
  NCORES     = c(30),
  batch_size = c(50000, 100000),
  {
    message(">>> Running NCORES=", NCORES, " batch_size=", batch_size)
    bench::mark(
      run_bench(NCORES, batch_size),
      iterations = 1,
      check = FALSE,
      memory = FALSE   # 👈 important for parallel code
    )
  }
)

print(results)

### 4kCpGs
#  expression   NCORES batch_size   min median `itr/sec` mem_alloc `gc/sec` n_itr
#1 run_bench(N…     30       5000 8.03m  8.03m   0.00207        NA     7.96     1
#2 run_bench(N…     30      10000 7.45m  7.45m   0.00224        NA     9.19     1

### 1000 CpGs:
#  expression   NCORES batch_size   min median `itr/sec` mem_alloc `gc/sec` n_itr
#1 run_bench(N…     20        500 2.75m  2.75m   0.00607        NA     4.21     1
#2 run_bench(N…     30        500  2.3m   2.3m   0.00725        NA     5.10     1
#3 run_bench(N…     20       1000 2.35m  2.35m   0.00708        NA     4.89     1
#4 run_bench(N…     30       1000 1.85m  1.85m   0.00902        NA     6.04     1

### 100 CpGs:
## A tibble: 9 × 15
#  expression  NCORES batch_size    min median `itr/sec` mem_alloc `gc/sec` n_itr
#  <bch:expr>   <dbl>      <dbl> <bch:> <bch:>     <dbl> <bch:byt>    <dbl> <int>
#1 run_bench(…      5         10  1.57m  1.57m    0.0106        NA    0.700     1
#2 run_bench(…     15         10  1.12m  1.12m    0.0148        NA    0.755     1
#3 run_bench(…     30         10  1.15m  1.15m    0.0145        NA    0.752     1
#4 run_bench(…      5         50  1.14m  1.14m    0.0147        NA    1.01      1
#5 run_bench(…     15         50 39.82s 39.82s    0.0251        NA    1.78      1
#6 run_bench(…     30         50 32.44s 32.44s    0.0308        NA    2.22      1
#7 run_bench(…      5        100  55.8s  55.8s    0.0179        NA    1.24      1
#8 run_bench(…     15        100 29.93s 29.93s    0.0334        NA    2.31      1
#9 run_bench(…     30        100 24.63s 24.63s    0.0406        NA    2.92      1

## Use <= 15G per thread
## A tibble: 9 × 15
#  expression   NCORES batch_size   min median `itr/sec` mem_alloc `gc/sec` n_itr
#  <bch:expr>    <dbl>      <dbl> <bch> <bch:>     <dbl> <bch:byt>    <dbl> <int>
#1 run_bench(N…     10        100 36.7s  36.7s    0.0272        NA     1.88     1
#2 run_bench(N…     20        100 28.4s  28.4s    0.0352        NA     2.46     1
#3 run_bench(N…     30        100 24.1s  24.1s    0.0415        NA     2.98     1
#4 run_bench(N…     10        500 37.1s  37.1s    0.0270        NA     1.86     1
#5 run_bench(N…     20        500 27.1s  27.1s    0.0368        NA     2.62     1
#6 run_bench(N…     30        500 24.4s  24.4s    0.0411        NA     3.00     1
#7 run_bench(N…     10       1000 37.6s  37.6s    0.0266        NA     1.83     1
#8 run_bench(N…     20       1000 27.5s  27.5s    0.0364        NA     2.62     1
#9 run_bench(N…     30       1000 25.2s  25.2s    0.0396        NA     2.97     1

##################################
## Comparisons WORKS (rm since) ##
##################################
##load("05_hvCpGalgorithm/resultsDir/Mariads/test/results_Maria_6906CpGs_0_8p0_0_65p1.RData")
##resV4 = results_Maria_6906CpGs_0_8p0_0_65p1
##load("05_hvCpGalgorithm/resultsDir/Mariads/results_Maria_pos2check_0_8p0_0_65p1.RData")
##resV3 = results_Maria_pos2check_0_8p0_0_65p1
##table(round(resV3,3) == round(resV4,3))
#### TRUE 
## 6906 

# ## Load hvCpG of Maria and matching mQTL
# cistrans_GoDMC_hvCpG_matched_control <- read.table("03_prepDatasetsMaria/cistrans_GoDMC_hvCpG_matched_control.txt", header = T)
# 
# ## cpg_names_all in here:
# prepData(analysis="Mariasarrays")
# length(cpg_names_all); head(cpg_names_all)
# ## [1] 406036
# # [1] "cg00000029" "cg00000108" "cg00000109" "cg00000165" "cg00000236" "cg00000289"
# 
# ## In which positions are they?
# cpg_names_all[cpg_names_all %in% cistrans_GoDMC_hvCpG_matched_control$hvCpG_name] %>% length
# cpg_names_all[cpg_names_all %in% cistrans_GoDMC_hvCpG_matched_control$controlCpG_name]  %>% length
# 
# sub_cistrans_GoDMC_hvCpG_matched_control <- cistrans_GoDMC_hvCpG_matched_control[
#   cistrans_GoDMC_hvCpG_matched_control$hvCpG_name %in% cpg_names_all &
#     cistrans_GoDMC_hvCpG_matched_control$controlCpG_name %in% cpg_names_all,]
# 
# ## Positions of targets in cpg_names:
# match(sub_cistrans_GoDMC_hvCpG_matched_control$hvCpG_name, cpg_names_all)
# match(sub_cistrans_GoDMC_hvCpG_matched_control$controlCpG_name, cpg_names_all)
# 
# ## test
# source_M_1CpG(cpgpos = match(sub_cistrans_GoDMC_hvCpG_matched_control$hvCpG_name, cpg_names_all)[2])
# 
# pos2check <- c(match(sub_cistrans_GoDMC_hvCpG_matched_control$hvCpG_name, cpg_names_all),
#                match(sub_cistrans_GoDMC_hvCpG_matched_control$controlCpG_name, cpg_names_all))

# ## Your grid
# params <- seq(0.4, 1, 0.05)
# grid <- expand.grid(p0 = params, p1 = params)
# 
# apply(grid, 1, function(row) {
#   runAndSave(analysis = "Maria", cpgPos_vec = pos2check, resultDir="05_hvCpGalgorithm/resultsDir/Mariads/",
#              NCORES=30, p0=as.numeric(row["p0"]), p1=as.numeric(row["p1"]))
# })

######################################################################
## Part 2. Subset Maria's data following Atlas data, to check power ##
######################################################################

## In Atlas data, a lot of CpG have only 3 datasets with 3 samples

###############
## Done before: python3 S02.prepare_REDUCED_arrays_mimicAtlas.py 
###############

for (s in 3:5){
    for (d in 3:30){
        if (file.exists(paste0("/home/alice/arraysh5_reducedMimicAtlas_", s, "samples_", d, "datasets/all_scaled_matrix.h5")))
        system.time(runAndSave(analysis = paste0("MariasarraysREDUCED_", s, "samples_", d, "datasets"), cpgPos_vec = pos2check,
                               resultDir="05_hvCpGalgorithm/resultsDir/Mariads/", NCORES=30, p0=0.80, p1=0.65))
    }
}

for (d in 3:30){
    if (file.exists(paste0("/home/alice/arraysh5_reducedMimicAtlas_allsamples_", d, "datasets/all_scaled_matrix.h5")))
        system.time(runAndSave(analysis = paste0("MariasarraysREDUCED_allsamples_", d, "datasets"), cpgPos_vec = pos2check,
                               resultDir="05_hvCpGalgorithm/resultsDir/Mariads/", NCORES=30, p0=0.80, p1=0.65))
}


###############
## To do after:
###############

## rm -rf /home/alice/arraysh5_reducedMimicAtlas_*samples*
