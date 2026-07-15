## Load hyperVarMeth (installed in my CS HPC home R.4.4.1 directory)                                                                                                                            
library(hyperVarMeth)                                                                                                                                                                         

################################################################################
## Run the algorithm on: full + reduced version of the array data where we kept 
## 3 people randomly chosen per dataset + reduced with 2 people
################################################################################

## Data must be prepared ahead (in 01_dataPrep)

## Output directory
result_dir <- "/home/alice/2024_hvCpG/B_MultiTissues/resultsDir_gitIgnored/Arrays"
message(paste0("If new, results will be saved in dir: ", result_dir))

runFull <- function(data_dir, cpg_names_all, analysis, p0, p1, minind=3){
  message(paste0("Run algo: ", analysis, " p0 = ", p0, "p1 = ", p1))
  
  system.time(hyperVarMeth::runAndSave_fast(
    analysis = analysis,
    cpg_names_vec = cpg_names_all,
    dataDir = data_dir,
    resultDir = result_dir,
    NCORES = 30,
    p0 = p0,
    p1 = p1,
    batch_size = 10000,
    minind = minind)
  )
}

## Run
runFull(data_dir = "/home/alice/arraysh5/",
        cpg_names_all = rhdf5::h5read("/home/alice/arraysh5/all_matrix_noscale.h5", "cpg_names"),
        analysis = "Arrays_all",
        p0=0.80, p1 = 0.65)

runFull(data_dir = "/home/alice/arraysh5_3ind/",
        cpg_names_all = rhdf5::h5read("/home/alice/arraysh5_3ind/all_matrix_noscale.h5", "cpg_names"),
        analysis = "Arrays_3ind",
        p0=0.80, p1 = 0.65)

runFull(data_dir = "/home/alice/arraysh5_2ind/",
        cpg_names_all = rhdf5::h5read("/home/alice/arraysh5_2ind/all_matrix_noscale.h5", "cpg_names"),
        analysis = "Arrays_2ind",
        p0=0.80, p1 = 0.65, minind = 2)

## With stricter p0/p1
runFull(data_dir = "/home/alice/arraysh5/",
        cpg_names_all = rhdf5::h5read("/home/alice/arraysh5/all_matrix_noscale.h5", "cpg_names"),
        analysis = "Arrays_all",
        p0=0.80, p1 = 0.9)

runFull(data_dir = "/home/alice/arraysh5_3ind/",
        cpg_names_all = rhdf5::h5read("/home/alice/arraysh5_3ind/all_matrix_noscale.h5", "cpg_names"),
        analysis = "Arrays_3ind",
        p0=0.80, p1 = 0.9)

runFull(data_dir = "/home/alice/arraysh5_2ind/",
        cpg_names_all = rhdf5::h5read("/home/alice/arraysh5_2ind/all_matrix_noscale.h5", "cpg_names"),
        analysis = "Arrays_2ind",
        p0=0.80, p1 = 0.9, minind = 2)
