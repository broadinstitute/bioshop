Bioinformatic prototyping library, currently focused on joint variant filtering.  Currently supports three commands, `etl`, `fit`, and `call`.
Use with a joint callset VCF by training a model with a high confidence allele set (such as dbSNP).  Example program flow:

1. Combine a joint callset VCF with a high confidence VCF, and use the ETL command to generate a dataframe of alleles.
2. Train a model with the `fit` command.
3. Use the trained model to `call` joint VCFs.

```
usage: newt etl [-h] --query_vcf QUERY_VCF_PATH --target_vcf TARGET_VCF_PATH
                [--assembly ASSEMBLY_NAME] [-o OUTPUT_PATH] -R REGION [-S STRAT_INTERVALS]

Train allele specific classification

options:
  -h, --help            show this help message and exit
  --query_vcf QUERY_VCF_PATH
                        Path to VCF to call
  --target_vcf TARGET_VCF_PATH
                        Path to VCF with valid calls from population
  --assembly ASSEMBLY_NAME
                        Name of the geome assembly to use
  -o OUTPUT_PATH, --output OUTPUT_PATH
                        Path for generated Pandas dataframe
  -R REGION, --region REGION
                        Region to generate results from
  -S STRAT_INTERVALS, --stratification STRAT_INTERVALS
                        Interval file for labeling lookup
```
  
```
  usage: newt fit [-h] -i INPUT_LIST [-o CLASSIFIER_PATH] [--classifier {rf,gb,mlp,kn,ada,xgb}]
                [--test-frac TEST_FRAC] [--random-seed RANDOM_SEED] [--combine]

Fit classifier to data

options:
  -h, --help            show this help message and exit
  -i INPUT_LIST, --input INPUT_LIST
                        Path to one or more panda dataframes
  -o CLASSIFIER_PATH, --output CLASSIFIER_PATH
                        Path to write classifier state
  --classifier {rf,gb,mlp,kn,ada,xgb}
                        Classifier to use
  --test-frac TEST_FRAC
                        Fraction to hold in reserve for testing
  --random-seed RANDOM_SEED
                        Random seed to use (randomly assigned if left unset)
  --combine             Generate two different models for snp/non-snp (default) or a single,
                        combined model
```
 
```
 usage: newt call [-h] --query_vcf QUERY_VCF_PATH [--classifier CLASSIFIER_PATH]
                 [--assembly ASSEMBLY_NAME] [-o OUTPUT_VCF_PATH] -R REGION [-S STRAT_INTERVALS]

Call VCF sites based on allelic training data

options:
  -h, --help            show this help message and exit
  --query_vcf QUERY_VCF_PATH
                        Path to VCF to call
  --classifier CLASSIFIER_PATH
                        Path to pickled classifier
  --assembly ASSEMBLY_NAME
                        Name of the geome assembly to use
  -o OUTPUT_VCF_PATH, --output OUTPUT_VCF_PATH
                        Path to generated VCF
  -R REGION, --region REGION
                        Region to generate results from
  -S STRAT_INTERVALS, --stratification STRAT_INTERVALS
                        Interval file for labeling lookup
```
