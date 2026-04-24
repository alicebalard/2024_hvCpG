## Run the algorithm on a reduced version of the array data, where we kept only 3 people randomly chosen per dataset

## Data directory
data_dir <- "/home/alice/arraysh5_3ind/"

## Output directory
result_dir <- "/home/alice/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Arrays"
message(paste0("If new, results will be saved in dir: ", result_dir))

## Load algorithm
source("/home/alice/2024_hvCpG/05_hvCpGalgorithm/hvCpG_algorithm_detection_v6.R")

## Names of all CpG in the dataset:                                                                                                                                                   
cpg_names_all <- rhdf5::h5read("/home/alice/arraysh5_3ind/all_matrix_noscale.h5", "cpg_names")

## Run
message("Run algo:")

system.time(runAndSave_fast(
    analysis = "Arrays_3indperds",
    cpg_names_vec = cpg_names_all,
    dataDir = data_dir,
    resultDir = result_dir,
    NCORES = 30,
    p0 = 0.80,
    p1 = 0.65,
    batch_size = 10000)
)
