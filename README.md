# Bioshop

bioshop is a bioinformatic prototyping library, currently focused on joint variant filtering.  Currently supports three commands, `etl`, `fit`, and `call`.
Use with a joint callset VCF by training a model with a high confidence allele set (such as dbSNP).  Example program flow:

1. Combine a joint callset VCF with a high confidence VCF, and use the `etl` command to generate a dataframe of alleles.
2. Train a model with the `fit` command.
3. Use the trained model to `call` joint VCFs.

### Installation

You can use `pip install` for docker for installation.

### Extract, translate and load the training data

If our VCF file was called `my_joint.vcf`, we had four intervals regions we want to 
label, and we want to use dbSNP as the allelic phonebook for all of chromosome 10, an example invocation would look like:

```
$ newt etl \
    --query_vcf /data/my_joint.vcf \
    --target_vcf /data/GCF_000001405.39.gz \
    -S homopol=/data/homopol.interval_list \
    -S dup=/data/duplicates.interval_list \
    -S high_gc=/data/high_gc.interval_list \
    -S low_complex=/data/low_complex.interval_list \
    -R chr10 \
    --output /data/my_joint_dataframe.pickle
 ```
 
 ### Fitting the model
 
 After this, we can generate our model.  This is an example of training a combined model (ie. SNPs and INDELs in one model) with xgboost:
 
 ```
$ newt fit \
      --combine \
      --classifier xgb \
      --input /data/my_joint_dataframe.pickle \
      --output /data/my_joint_model.pickle
```

### Calling VCFs

Once we generate our model, we can now call new joint VCFs.
 
```
$ newt call \
    --query_vcf /data/my_other_joint.vcf \
    --classifier /data/my_joint_model.pickle \
    -S homopol=/data/homopol.interval_list \
    -S dup=/data/duplicates.interval_list \
    -S high_gc=/data/high_gc.interval_list \
    -S low_complex=/data/low_complex.interval_list \
    -R chr1 \
    --output /data/my_other_joint_called.vcf
```

### --help for help

Use `--help` for command line help.  For example:

```
$ newt etl --help
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
