#nNIL introgression SNP comparisons on NCSU HPC

##change directory to hpc storage
cd /rs1/researchers/j/jholland/nNILs

##install a conda environment with samtools and bcftools
#conda create --name=snp_filtering samtools bcftools
conda activate snp_filtering

## i tried ncsu spack bioinfo modules, but doesn't have what i need
#module use --append /usr/local/usrapps/bioinfo/modulefiles
#module avail


##I downloaded SNP calls from Grzybowski maize Plant Journal 2023
#cd  Grzybowski/

##get sample names in Grzybowski
#bcftools query -l chr_1_imputed.vcf.gz > grzybowski_samples.txt

##Turns out that CML103 and Mo18W are MISSING from Grzybowski data set.

##Originally, I made a comprehensive list of nNIL SNP positions in V5 coordinates, but see below, we actually want them all in V3 to compare to HapMap3 data
#combine the gbs and chip SNP positions on nNILs into a single bed file
#note some positions will be included twice, I don't think that is a problem for what we are going to do
#cat nNIL_chip_SNP_positions_converted_to_V5.bed nNIL_gbs_SNP_positions_converted_to_V5.bed > nNIL_chip_gbs_combined_SNP_positions_converted_to_V5.bed


##So, instead, will use the HapMap 3 data from Bukowski et al 2017 paper
##Which are in version AGP 3 coordinates!!
##Cyverse path to these data is:
##/iplant/home/shared/commons_repo/curated/Qi_Sun_Zea_mays_haplotype_map_2018/282_onHmp321


##Instead of uplifting SNP positions to V5, used Ensembl to convert Laura's GBS positions from V4 to V3 !
##And subset the HapMap3 SNPs to the overlapping set of positions and only for the nNIL founders

#Get the list of HapMap3 samples
#cd HapMap3
#bcftools query -l c1_282_corrected_onHmp321.vcf.gz > hapmap3_samples.txt


#I generated a list of the nNIL GBS SNPs projected to AGP v3 coordinates to compare to HapMap3
# See R script "SNP_info_nNILs.Rmd"
#subset each chromosome vcf to the subset of matching GBS SNPs and founders (we will do chip SNPs separately)
#require compressing (even though they say they are gz they are not) and indexing first
#loop over chromosomes
for i in {1..10}
    do
       #bcftools view -Oz -o "c${i}_282_corrected_bgz_onHmp321.vcf.gz" "c${i}_282_corrected_onHmp321.vcf.gz"
       #tabix "c${i}_282_corrected_bgz_onHmp321.vcf.gz"
       bcftools view -S nNIL_founder_samples.txt -R nNIL_gbs_SNPs_v3.bed -Oz "c${i}_282_corrected_bgz_onHmp321.vcf.gz" > "c${i}_nNIL_founders_gbsSNPs_Hmp321.vcf.gz"
    done

#get the HapMap 3 AGP_V3 SNP positions for each chromosome and concat into one file (use >> to append)
for i in {1..10}
    do
      bcftools query -f '%CHROM %POS %REF %ALT\n' "c${i}_282_corrected_bgz_onHmp321.vcf.gz" >> SNPs_282_Hmp321.txt
    done

# combine all ten chromosomes, keep only biallelic SNPs, and check how many SNPs are retained in each subset chr vcf
bcftools concat c1_nNIL_founders_gbsSNPs_Hmp321.vcf.gz c2_nNIL_founders_gbsSNPs_Hmp321.vcf.gz \
c3_nNIL_founders_gbsSNPs_Hmp321.vcf.gz c4_nNIL_founders_gbsSNPs_Hmp321.vcf.gz c5_nNIL_founders_gbsSNPs_Hmp321.vcf.gz c6_nNIL_founders_gbsSNPs_Hmp321.vcf.gz \
c7_nNIL_founders_gbsSNPs_Hmp321.vcf.gz c8_nNIL_founders_gbsSNPs_Hmp321.vcf.gz c9_nNIL_founders_gbsSNPs_Hmp321.vcf.gz c10_nNIL_founders_gbsSNPs_Hmp321.vcf.gz | \
bcftools view --types snps -m 2 -M 2 -Oz -o nNIL_founders_gbsSNPs_Hmp321_AGPv3.vcf.gz

bcftools stats nNIL_founders_gbsSNPs_Hmp321_AGPv3.vcf.gz > nNIL_founders_gbsSNPs_Hmp321_AGPv3.stats.txt

#output the info for the remaining SNPs, compare reference/major/minor to the GBS set in R
bcftools query -f '%CHROM %POS %REF %ALT\n' nNIL_founders_gbsSNPs_Hmp321_AGPv3.vcf.gz > nNIL_founder_Hmp321SNPs_overlapping_GBS.txt

#get the info on the GBS SNPs
bcftools view -Oz -o SNPs_orig_gbs_V4_nNILs.vcf.gz SNPs_orig_gbs_V4_nNILs.vcf
tabix SNPs_orig_gbs_V4_nNILs.vcf.gz
bcftools query -f '%CHROM %POS %REF %ALT\n' SNPs_orig_gbs_V4_nNILs.vcf.gz > SNPs_GBS_nNIL_v4.txt

#Used R script "SNP_info_nNILs.Rmd" to generate a final bed file that contains only the consistent subset of SNPs
#also note I forgot to include -S nNIL_founder_samples.txt above, so have to do it here...
#these are SNPs where HapMap REF allele = gbs major allele and HapMap ALT allele = gbs minor allele
tabix nNIL_founders_gbsSNPs_Hmp321_AGPv3.vcf.gz
bcftools view -S nNIL_founder_samples.txt -R nNIL_gbs_SNPs_congruent_w_HapMap_v3.bed nNIL_founders_gbsSNPs_Hmp321_AGPv3.vcf.gz  |\
  bcftools sort -Oz -o nNIL_founders_consistent_gbsSNPs_v3.vcf.gz

#check it
bcftools stats nNIL_founders_consistent_gbsSNPs_v3.vcf.gz > nNIL_founders_consistent_gbsSNPs_v3.stats.txt

#now subset the original GBS vcf of nNILs to the same set of consistent markers
cd ..

###WAIT! FIRST CONVERT GBS VCF TO V3 POSITIONS!!!!

bcftools view -Oz -o SNPs_orig_gbs_V4_nNILs.vcf.gz SNPs_orig_gbs_V4_nNILs.vcf
tabix SNPs_orig_gbs_V4_nNILs.vcf.gz
bcftools view -R nNIL_gbs_SNPs_congruent_w_HapMap_v3.bed -Oz -o orig .... gbs_V4_nNILs.vcf.gz orig_gbs_V4_nNILs.vcf.gz

#############################################################
# Make subset of HapMap GBS file for nNIL founders
# and markers used for chip genotyping
# to compare chip genotypes and hapmap genotypes
#############################################################


for i in {1..10}
    do
       #bcftools view -Oz -o "c${i}_282_corrected_bgz_onHmp321.vcf.gz" "c${i}_282_corrected_onHmp321.vcf.gz"
       #tabix "c${i}_282_corrected_bgz_onHmp321.vcf.gz"
       bcftools view -S nNIL_founder_samples.txt -R nNIL_chip_SNPs_v3.bed -Oz "c${i}_282_corrected_bgz_onHmp321.vcf.gz" > "c${i}_nNIL_founders_chipSNPs_Hmp321.vcf.gz"
    done

# combine all ten chromosomes, keep only biallelic SNPs, and check how many SNPs are retained in each subset chr vcf
bcftools concat c1_nNIL_founders_chipSNPs_Hmp321.vcf.gz c2_nNIL_founders_chipSNPs_Hmp321.vcf.gz \
c3_nNIL_founders_chipSNPs_Hmp321.vcf.gz c4_nNIL_founders_chipSNPs_Hmp321.vcf.gz c5_nNIL_founders_chipSNPs_Hmp321.vcf.gz c6_nNIL_founders_chipSNPs_Hmp321.vcf.gz \
c7_nNIL_founders_chipSNPs_Hmp321.vcf.gz c8_nNIL_founders_chipSNPs_Hmp321.vcf.gz c9_nNIL_founders_chipSNPs_Hmp321.vcf.gz c10_nNIL_founders_chipSNPs_Hmp321.vcf.gz | \
bcftools view --types snps -m 2 -M 2 -Oz -o nNIL_founders_chipSNPs_Hmp321_AGPv3.vcf.gz

bcftools stats nNIL_founders_chipSNPs_Hmp321_AGPv3.vcf.gz > nNIL_founders_chipSNPs_Hmp321_AGPv3.stats.txt

#############################################################
# Finally, go back to original HapMap files and 
# filter to the nNIL founders and sites with <= 20% missing data
# And zero hets among the founder subset
# These files will be transformed to numericalized genotypes in tassel gui
#############################################################

# I am assuming that missing data proportions need to be recomputed AFTER subsetting
# Not sure this is true, but just in case we will subset, fill tags, 
# then filter to include only biallelic snps with no missing data, no hets, and monomorphic
# -m2 -M2 gets the biallelic markers (min (m) and max (M) number of alleles = 2)
# -v snps gets the SNPs
# -c 1 gets polymorphic markers in the subset
# -g ^het gets markers with NO hets
# also fix the sample names using file nNIL_founder_samples_rename.txt

for i in {1..10}
    do 
    echo chrom $i
    bcftools view -S nNIL_founder_samples.txt -Ou "c${i}_282_corrected_bgz_onHmp321.vcf.gz" | \
    bcftools reheader -s nNIL_founder_samples_rename.txt | \
    bcftools +fill-tags -Ou -- -t F_MISSING | \
    bcftools view -m2 -M2 -v snps -c 1 -i 'F_MISSING<0.01' -g ^het -Oz -o "c${i}_nNIL_founders_miss0_Hmp321.vcf.gz"
    done  


# make a file with all the file names, to get the ordering correct we have to put c10 at end by hand
ls *nNIL_founders_miss0_Hmp321.vcf.gz | tail -n +2  > miss0_filenames.txt #this gets chr 10, 1 - 9, then drops first line
echo c10_nNIL_founders_miss0_Hmp321.vcf.gz >> miss0_filenames.txt #this adds chr 10 at bottom

# concatenate the subset vcfs into a common vcf
bcftools concat -f miss0_filenames.txt -Oz -o nNIL_founders_miss0_Hmp321.vcf.gz