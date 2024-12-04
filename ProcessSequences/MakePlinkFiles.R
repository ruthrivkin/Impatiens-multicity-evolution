
setwd("/Users/ruthrivkin/Dropbox/Grad School 2015-2020/Impatiens capensis/Data files/9CitiesSeq/Population Genetics/NewAssembly/")


#Create plink files
system("~/plink/plink --vcf variants_hardfilter.dp10.geno.maf.vcf --recode --allow-extra-chr --out variants_hardfilter.dp10.geno.maf")
#Total genotyping rate is 0.963188.
#9021 variants and 228 people pass filters and QC.
  
#Update ids to city and population level
system("~/plink/plink --file variants_hardfilter.dp10.geno.maf --allow-extra-chr --update-ids UpdateID_city.txt --recode --out variants_hardfilter.dp10.geno.maf_city")
system("~/plink/plink --file variants_hardfilter.dp10.geno.maf --allow-extra-chr --update-ids UpdateID_pop.txt --recode --out variants_hardfilter.dp10.geno.maf_pop")

#Calculate ld - city only for now
#update map first
system("awk '{$2=FNR-1; print}' OFS='\t' variants_hardfilter.dp10.geno.maf_city.map > updated.loci.map")
system("~/plink/plink --ped variants_hardfilter.dp10.geno.maf_city.ped --map updated.loci.map --allow-extra-chr --indep-pairwise 10 1 0.1 --out snps_ld") #5056 variants removed

#Remove snps in ld
system("~/plink/plink --ped variants_hardfilter.dp10.geno.maf_city.ped --map updated.loci.map --allow-extra-chr --exclude snps_ld.prune.out --recode --out variants_hardfilter.dp10.geno.maf.ld_city")
#3965 variants and 228 people pass filters and QC

#Lastly filter for hwe
system("~/plink/plink --file variants_hardfilter.dp10.geno.maf.ld_city --allow-extra-chr --hwe 0.001 --recode --out variants_hardfilter.dp10.geno.maf.ld.hwe_city")
#3261 variants and 228 people pass filters and QC.


#Convert files to raw format
#hwe and ld
system("~/plink/plink --file variants_hardfilter.dp10.geno.maf.ld_city --allow-extra-chr --hwe 0.001 --recode A --out variants_hardfilter.dp10.geno.maf.ld.hwe_city")

#LD only
system("~/plink/plink --ped variants_hardfilter.dp10.geno.maf_city.ped --map updated.loci.map --allow-extra-chr --exclude snps_ld.prune.out --recode A --out variants_hardfilter.dp10.geno.maf.ld_city")

#Also make raw files for population level analyses (fst)
system("~/plink/plink --file variants_hardfilter.dp10.geno.maf.ld_city -update-ids UpdateID_pop.txt --allow-extra-chr --hwe 0.001 --recode A --out variants_hardfilter.dp10.geno.maf.ld.hwe_pop")

