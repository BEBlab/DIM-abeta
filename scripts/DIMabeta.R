
DIMabeta<-function(wkd, dimsum_fitness_replicates_file, dimsum_variant_data_merge_file){


	# set working directory, which contains a folder with required data and a folder with all R scripts
	setwd(wkd)

	#load DiMSum output data
	load(dimsum_fitness_replicates_file)
	load(dimsum_variant_data_merge_file)


	source("scripts/01_processed data file.R")
	source("scripts/02_FDR and mutation types.R")
	source("scripts/03_DMS data quality.R")
	source("scripts/04_NS distributions.R")
	source("scripts/05_ROC.R")
	source("scripts/06_single aa mutations.R")
	source("scripts/07_dendogram.R")
	source("scripts/08_mutation types comparison.R")
	source("scripts/09_predictors.R")
	source("scripts/10_deletions.R")
	source("scripts/11_alternative cores.R")
	source("scripts/12_alternative Nterminus.R")
	source("scripts/13_truncations.R")
	source("scripts/14_FDR heatmaps.R")
	source("scripts/15_top nucleators.R")
	source("scripts/16_chimera.R")


}