## Run the algorithm on a reduced version of the array data, where we kept only 2 people randomly chosen per dataset

## Load hyperVarMeth (installed in my CS HPC home R.4.2 directory)                                                                                                                            
library(hyperVarMeth)                                                                                                                                                                         

## Data directory
data_dir <- "/home/alice/arraysh5_2ind/"

## Output directory
result_dir <- "/home/alice/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Arrays"
message(paste0("If new, results will be saved in dir: ", result_dir))

## Names of all CpG in the dataset:                                                                                                                                                   
cpg_names_all <- rhdf5::h5read("/home/alice/arraysh5_2ind/all_matrix_noscale.h5", "cpg_names")

## Run
message("Run algo:")

system.time(hyperVarMeth::runAndSave_fast(
    analysis = "Arrays_2indperds",
    cpg_names_vec = cpg_names_all,
    dataDir = data_dir,
    resultDir = result_dir,
    NCORES = 30,
    p0 = 0.80,
    p1 = 0.65,
    batch_size = 10000, minind = 2) ### Important!
    )
