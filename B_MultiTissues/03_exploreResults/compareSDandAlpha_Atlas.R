library(here)
library(dplyr)

## Compare standard deviation and alpha (calculated in paleo project)
fullres_Atlas10X <- readRDS(here("gitignore/fullres_Atlas10X"))
cpg_std_values <- read.csv(here("gitignore/cpg_std_values.csv"))
cpg_std_values <- cpg_std_values %>% dplyr::rename(name = CpG_Name)

head(fullres_Atlas10X)
head(cpg_std_values)

df <- merge(fullres_Atlas10X, cpg_std_values)
head(df)

rm(fullres_Atlas10X, cpg_std_values)

library(ggplot2)

dfplot <- df[sample(nrow(df), 100000), ]

ggplot(df, aes(x = alpha, y =Standard_Deviation)) +
  geom_point(data = dfplot, alpha = .3) +
  geom_smooth() 
