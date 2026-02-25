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
```



