## Run the algorithm on either the full or a reduced version of the array data,
## where we kept only N people randomly chosen per dataset

## Load hyperVarMeth (installed in my CS HPC home R.4.4.1 directory)
library(hyperVarMeth)

## Output directory
result_dir <- "/home/alice/2024_hvCpG/B_MultiTissues/resultsDir_gitIgnored"
message(paste0("If new, results will be saved in dir: ", result_dir))

runAll <- function(data_dir, analysis, n = "all", p0, p1){

    if (!exist(data_dir)){
        print(paste0("No ", data_dir, " found"))
    } else {        
        ## Names of all CpG in the dataset:
        cpg_names_all <- rhdf5::h5read(paste0(data_dir, "all_matrix_noscale.h5"), "cpg_names")

        message("Run algo:")
        if (n == 2){
            minind = 2
        } else {
            minind = 3
        }
        
        system.time(hyperVarMeth::runAndSave_fast(
                                      analysis = analysis,
                                      cpg_names_vec = cpg_names_all,
                                      dataDir = data_dir,
                                      resultDir = result_dir,
                                      NCORES = 30,
                                      p0 = p0,
                                      p1 = 01,
                                      batch_size = 10000, minind = minind)
                    )
    }
}

data_dir = paste0("/home/alice/arraysh5_all")
runAll(data_dir, n = i, paste0("Arrays_all"), p0=0.8, p1=0.65)

for (i in c(2, 3)){
    data_dir = paste0("/home/alice/arraysh5_", i, "ind/")
    runAll(data_dir, n = i, paste0("Arrays_", i, "indperds"), p0=0.8, p1=0.65)
}

## Stricter parameters
data_dir = paste0("/home/alice/arraysh5_all")
runAll(data_dir, n = i, paste0("Arrays_all"), p0=0.8, p1=0.9)

for (i in c(2, 3, 4, 5, 6, 10)){
    runAll(data_dir, n = i, p0=0.8, p1=0.9)
}
