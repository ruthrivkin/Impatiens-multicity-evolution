#!/usr/bin/env Rscript
#SBATCH -J pixyplot
#SBATCH -o logs/pixyp.out
#SBATCH -e logs/pixyp.err
#SBATCH -n 1
#SBATCH -t 10000
#SBATCH --mem=200G
#SBATCH --partition=standard
#SBATCH --account=grdi_genarcc
#SBATCH --comment="image=registry.maze.science.gc.ca/ssc-hpcs/generic-job:ubuntu22.04,tmpfs_size=20G"

.libPaths(new = "/gpfs/fs7/grdi/genarcc/wp3/PolarBears/env/R")
library(cli)
library(rlang)

setwd("/gpfs/fs7/grdi/genarcc/wp3/PolarBears/ICurb/Analysis/Pixy")
getwd()

library(tidyverse)

pixy_to_long <- function(pixy_files){

  pixy_df <- list()

  for(i in 1:length(pixy_files)){

    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])

    if(stat_file_type == "pi"){

      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)

      pixy_df[[i]] <- df


    } else{

      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value")
      pixy_df[[i]] <- df

    }

  }

  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)

}

pixy_folder <- "Output/pop"
pixy_files <- list.files(pixy_folder, full.names = TRUE)
pixy_df <- pixy_to_long(pixy_files)

pixy_df_scaffold <- pixy_df %>%
  separate(chromosome, c("sc", "scaffold", "ex"),  "_", extra = "merge") %>%
  select(!c(sc, ex))

saveRDS(pixy_df_scaffold, "pixy_scaffold.rds")