# Impatiens-multicity-evolution
Code for analysis of neutral and adaptive evolution of Impatiens capensis in 10 cities in Ontario

There are six folders that group analyses together: 
1. ProcessSequences: Trim, assemble, and call SNPs across samples. Also create plink files for downstream analyses and modified GBS-preprocess files from https://github.com/relshire/GBS-PreProcess
2. GeneticDiversity: Calculate genetic diversity statistics and run linear models of habitat change on genetic diversity
3. GeneticStructure: Estimate population structure at broad and local scales
4. DemographicHistory: Calculate current Ne and and recent changes in demographic species
5. GenotypebyEnv: Run GEA models associating habitat with putatively adaptive loci
6. Scripts: bash scripts formatted to run on a slurm workload manager on an HPC. These files link to different files within the other folders.

Note folder organization in Github is not the same as in the code because the work was split across an HPC and personal computer with differing file structures. 
