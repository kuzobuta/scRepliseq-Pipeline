#!/bin/bash
#Setup the path for local run
#This command do not need to run for docker

local_path=`pwd`

#file set1 under scripts
files=("Step3_load_mapped_reads.sh" \
"Step4_MAD_score_QC.sh" \
"Step5_1_Check_G1_cells.sh" \
"Step5_2_Merge_G1_cells.sh" \
"Step6_log2_median_RT_scores.sh" \
"Step7_Binarization.sh")

for file in ${files[@]}
do
(rm "scripts/$file" && sed -e "s/\/usr\/local\/bin/\$local_path/g" | awk -v d="${local_path}" '{if(NR==1){print $0"\nlocal_path=""\""d"\""}else{print}}' > scripts/$file) < "scripts/$file"
done

#file set2 under util
files=("Fragment2bin.sh" \
"Predict_SEQXE_Primer.sh")

for file in ${files[@]}
do
(rm "util/$file" && sed -e "s/\/usr\/local\/bin/\$local_path/g" | awk -v d="${local_path}" '{if(NR==1){print $0"\nlocal_path=""\""d"\""}else{print}}' > util/$file) < "util/$file"
done

#file set3 (Rscripts)
files=("Generate_IGV_session_two_data.R" \
"Generate_IGV_session_one_data.R" \
"Step6_R_log2_median_RT_scores_HMM.R" \
"Step6_R_log2_median_RT_scores.R" \
"Step7_R_Binarization.R")

for file in ${files[@]}
do
(rm "util/$file" && awk -v d="${local_path}" '{if(NR==2){print "path=""\""d"\""}else{print}}' > util/$file) < "util/$file"
done
