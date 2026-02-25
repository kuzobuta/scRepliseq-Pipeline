# scRepliseq-Pipeline version 260225

Release notes can be found [here](https://github.com/kuzobuta/scRepliseq-Pipeline/blob/master/Change_logs.md).

Citation [ Miura et al., Nature Protocols, 2020 ](https://www.nature.com/articles/s41596-020-0378-5)

## Test run
To perform a test run of this pipeline, below is the commands:

```
## Install the conda environment
conda create -f sample_data/env/screpliseq_v1.5.yml
# or
mamba env create -f sample_data/env/screpliseq_v1.5.yml

## Run the test scripts 
conda activate screpliseq_v1.5
bash test_run.sh &> test.run.log

## test.run.log
# Start test run:  Wed Feb 25 11:52:29 JST 2026
# Done step0:      Wed Feb 25 11:52:30 JST 2026
# Done FastQC:     Wed Feb 25 11:52:34 JST 2026
# Done step1:      Wed Feb 25 11:52:38 JST 2026
# Done step2:      Wed Feb 25 11:52:56 JST 2026
# Done step3:      Wed Feb 25 11:55:37 JST 2026
# Done step4:      Wed Feb 25 11:55:46 JST 2026
# Done step5-1:    Wed Feb 25 11:57:59 JST 2026
# Done step5-2:    Wed Feb 25 11:58:17 JST 2026
# Done tagDensity: Wed Feb 25 12:00:54 JST 2026
# Done step6:      Wed Feb 25 12:02:02 JST 2026
# Done step7:      Wed Feb 25 12:04:34 JST 2026
# End   test run:  Wed Feb 25 12:04:34 JST 2026
# Test run finished successfully.
##
## Elapsed time: 12 min 04 sec
##
```



