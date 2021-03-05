# Population-based BrdU-IP analysis
---

**Before performing BrdUIP data analysis, please geenerate the bam file or fragment R file for each Early-S or Late-S data**

* To generate fragment R file, you may use our scRepli-seq pipeline from Step1 to Step3.

## Example script for fragment R file

* Using sample_data 

```
#Usage: StepX1_BrdUIP_from_R_file.sh [Early_S] [Late_S] [name] [out_dir] [blacklist] [genome_file]
#[Early_S] [Late_S] are fragment R file

./StepX1_BrdUIP_from_R_file.sh \
../sample_data/fragment_data_for_BrdUIP/1_RPE1_IP_EarlyS_mapq10_blacklist_fragment.Rdata \
../sample_data/fragment_data_for_BrdUIP/2_RPE1_IP_LateS_mapq10_blacklist_fragment.Rdata \
RPE1_IP_test \
test_results/ \
../sample_data/blacklist/hg19-blacklist-v1_id.bed \
../sample_data/genome_file/UCSC_hg19_female.fa.fai
```

* Optionally, from bam file

```
#Usage: StepX1_BrdUIP_from_bam_file.sh [Early_S] [Late_S] [name] [out_dir] [blacklist] [genome_file]
#[Early_S] [Late_S] are bam file

./StepX1_BrdUIP_from_bam_file.sh \
[Path to EarlyS bam file] \
[Path to LateS bam file] \
RPE1_IP_test \
test_results/ \
../sample_data/blacklist/hg19-blacklist-v1_id.bed \
../sample_data/genome_file/UCSC_hg19_female.fa.fai
```




You can see the processed data using `sample_data/fragment_data_for_BrdUIP/ ` under `test_results`.
