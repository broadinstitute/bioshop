version 1.0

# This is a workflow for filtering a joint callset VCF using INFO level annotations (so filtering is at the site level).
# Note that the input VCFs here may be sharded by genomic position which may be helpful for large cohorts. The script
# will output the same number of shards that are input.
# This portion of the filtering pipeline will assign a SCORE INFO field annotation to each site, but does not yet apply
# the filtering threshold to the final VCF.

workflow JointVcfFiltering {
  input {
    File input_vcf
    File input_vcf_index
    File sites_only_vcf
    File sites_only_vcf_index
    File unpadded_intervals
    String basename

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    String model_backend
    File? training_python_script
    File? scoring_python_script
    File? hyperparameters_json

    String gatk_docker
    File? extract_interval_list
    File? score_interval_list

    String snp_annotations
    String indel_annotations
    File? gatk_override

    Boolean use_allele_specific_annotations

    String snp_resource_args = "--resource:hapmap,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz --resource:omni,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz --resource:1000G,training=true,calibration=false gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    String indel_resource_args = "--resource:mills,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

    Int scatter_count = 20
    Int small_disk = 128
    Int medium_disk = 256
    Int large_disk = 1024
    Int huge_disk = 2048
  }

  parameter_meta {
    vcf: "VCF to call"
    sites_only_vcf: "VCF for training"
    basename: "Desired output file basename."
  }

  call ExtractVariantAnnotations as ExtractVariantAnnotationsSNPs {
    input:
      input_vcf = sites_only_vcf,
      input_vcf_index = sites_only_vcf_index,
      mode = "SNP",
      annotations = snp_annotations,
      resource_args = snp_resource_args,
      basename = basename,
      interval_list = extract_interval_list,
      use_allele_specific_annotations = use_allele_specific_annotations,
      gatk_override = gatk_override,
      gatk_docker = gatk_docker
  }

  call ExtractVariantAnnotations as ExtractVariantAnnotationsINDELs {
    input:
      input_vcf = sites_only_vcf,
      input_vcf_index = sites_only_vcf_index,
      mode = "INDEL",
      annotations = indel_annotations,
      resource_args = indel_resource_args,
      basename = basename,
      interval_list = extract_interval_list,
      use_allele_specific_annotations = use_allele_specific_annotations,
      gatk_override = gatk_override,
      gatk_docker = gatk_docker
  }

  call TrainVariantAnnotationModel as TrainVariantAnnotationModelSNPs {
    input:
      annots = ExtractVariantAnnotationsSNPs.annots,
      basename = basename,
      mode = "snp",
      model_backend = model_backend,
      python_script = training_python_script,
      hyperparameters_json = hyperparameters_json,
      gatk_override = gatk_override,
      gatk_docker = gatk_docker
  }

  call TrainVariantAnnotationModel as TrainVariantAnnotationModelINDELs {
    input:
      annots = ExtractVariantAnnotationsINDELs.annots,
      basename = basename,
      mode = "indel",
      model_backend = model_backend,
      python_script = training_python_script,
      hyperparameters_json = hyperparameters_json,
      gatk_override = gatk_override,
      gatk_docker = gatk_docker
  }

  call SplitIntervalList as interval_shards {
    input:
      interval_list = unpadded_intervals,
      scatter_count = scatter_count,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      disk_size = small_disk,
  }

  call SelectVariants as subset_input_vcf {
    input:
      input_intervals = interval_shards.output_intervals,
      input_vcf_file = input_vcf,
      input_vcf_index_file = input_vcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
  }

  scatter(idx in range(length(subset_input_vcf.output_vcf))) {
    call ScoreVariantAnnotations as ScoreVariantAnnotationsSNPs {
      input:
        vcf = subset_input_vcf.output_vcf[idx],
        vcf_index = subset_input_vcf.output_vcf_index[idx],
        basename = basename,
        mode = "SNP",
        model_backend = model_backend,
        python_script = scoring_python_script,
        annotations = snp_annotations,
        extracted_training_vcf = ExtractVariantAnnotationsSNPs.extracted_training_vcf,
        extracted_training_vcf_index = ExtractVariantAnnotationsSNPs.extracted_training_vcf_index,
        interval_list = score_interval_list,
        model_files = TrainVariantAnnotationModelSNPs.outputs,
        resource_args = snp_resource_args,
        use_allele_specific_annotations = use_allele_specific_annotations,
        gatk_override = gatk_override,
        gatk_docker = gatk_docker
    }

    call ScoreVariantAnnotations as ScoreVariantAnnotationsINDELs {
      input:
        vcf = ScoreVariantAnnotationsSNPs.output_vcf,
        vcf_index = ScoreVariantAnnotationsSNPs.output_vcf_index,
        basename = basename,
        mode = "INDEL",
        model_backend = model_backend,
        python_script = scoring_python_script,
        annotations = indel_annotations,
        extracted_training_vcf = ExtractVariantAnnotationsINDELs.extracted_training_vcf,
        extracted_training_vcf_index = ExtractVariantAnnotationsINDELs.extracted_training_vcf_index,
        interval_list = score_interval_list,
        model_files = TrainVariantAnnotationModelINDELs.outputs,
        resource_args = indel_resource_args,
        use_allele_specific_annotations = use_allele_specific_annotations,
        gatk_override = gatk_override,
        gatk_docker = gatk_docker
    }

  }

  output {
    Array[File] variant_scored_vcf = ScoreVariantAnnotationsINDELs.output_vcf
    Array[File] variant_scored_vcf_index = ScoreVariantAnnotationsINDELs.output_vcf_index
  }
}

task ExtractVariantAnnotations {
  input {
    String gatk_docker
    File? gatk_override
    File input_vcf
    File input_vcf_index
    String basename
    String mode
    String annotations
    String resource_args
    File? interval_list
    Boolean use_allele_specific_annotations

    Int memory_mb = 14000
    Int command_mem = memory_mb - 1000
  }
  Int disk_size = ceil(size(input_vcf, "GB") + 50)
  command {
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    gatk --java-options "-Xmx~{command_mem}m" \
      ExtractVariantAnnotations \
      -V ~{input_vcf} \
      -O ~{basename}.~{mode} \
      ~{annotations} \
      ~{if use_allele_specific_annotations then "--use-allele-specific-annotations" else ""} \
      ~{"-L " + interval_list} \
      --mode ~{mode} \
      ~{resource_args}

    if [ ! -f "~{basename}.~{mode}.vcf.gz.tbi" ]; then
      gatk --java-options "-Xmx~{command_mem}m" \
        IndexFeatureFile \
        -I "~{basename}.~{mode}.vcf.gz" \
        -O "~{basename}.~{mode}.vcf.gz.tbi"
    fi
  }
  output {
    File annots = "~{basename}.~{mode}.annot.hdf5"
    File extracted_training_vcf = "~{basename}.~{mode}.vcf.gz"
    File extracted_training_vcf_index = "~{basename}.~{mode}.vcf.gz.tbi"
    Array[File] outputs = glob("~{basename}.~{mode}.*")
  }
  runtime {
    docker: gatk_docker
    disks: "local-disk " + disk_size + " LOCAL"
    memory: memory_mb + " MiB"
  }
}

task TrainVariantAnnotationModel {
  input {
    String gatk_docker
    File? gatk_override
    File annots
    String basename
    String mode
    String model_backend
    File? python_script
    File? hyperparameters_json

    Int memory_mb = 14000
    Int command_mem = memory_mb - 1000
  }
  Int disk_size = ceil(size(annots, "GB") + 100)
  command <<<
    set -e

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    mode=$(echo "~{mode}" | awk '{print toupper($0)}')

    gatk --java-options "-Xmx~{command_mem}m" \
      TrainVariantAnnotationsModel \
      --annotations-hdf5 ~{annots} \
      -O ~{basename} \
      --model-backend ~{model_backend} \
      ~{"--python-script " + python_script} \
      ~{"--hyperparameters-json " + hyperparameters_json} \
      --mode $mode

  >>>
  output {
    Array[File] outputs = glob("~{basename}.~{mode}.*")
  }
  runtime {
    docker: gatk_docker
    disks: "local-disk " + disk_size + " LOCAL"
    memory: memory_mb + " MiB"
  }
}

task ScoreVariantAnnotations {
  input {
    String gatk_docker
    File? gatk_override
    File vcf
    File vcf_index
    String basename
    String mode
    String model_backend
    File? python_script
    String annotations
    String resource_args
    File extracted_training_vcf
    File extracted_training_vcf_index
    File? interval_list
    Array[File] model_files
    Boolean use_allele_specific_annotations

    Int memory_mb = 16000
    Int command_mem = memory_mb - 1000
  }
  Int disk_size = ceil(size(vcf, "GB") * 4 + 50)

  command {
    zgrep -v '#' ~{vcf} > empty.txt
    set -e

    if [ -s empty.txt ]; then
      ln -s ~{sep=" . && ln -s " model_files} .

      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

      gatk --java-options "-Xmx~{command_mem}m" \
        ScoreVariantAnnotations \
        ~{"-L " + interval_list} \
        -V ~{vcf} \
        -O ~{basename}.~{mode} \
        --model-backend ~{model_backend} \
        ~{"--python-script " + python_script} \
        --model-prefix ~{basename} \
        ~{annotations} \
        ~{if use_allele_specific_annotations then "--use-allele-specific-annotations" else ""} \
        -mode ~{mode} \
        --resource:extracted,extracted=true ~{extracted_training_vcf} \
        ~{resource_args}
    else
      echo "Input VCF was empty so we'll return the same VCF that was input."
      echo "Scores and annot hdf5 files will not be produced since the input was empty."
      ln -s ~{vcf} ~{basename}.~{mode}.vcf.gz
      ln -s ~{vcf_index} ~{basename}.~{mode}.vcf.gz.tbi
    fi

    if [ ! -f "~{basename}.~{mode}.vcf.gz.tbi" ]; then
      gatk --java-options "-Xmx~{command_mem}m" \
        IndexFeatureFile \
        -I "~{basename}.~{mode}.vcf.gz" \
        -O "~{basename}.~{mode}.vcf.gz.tbi"
    fi
  }
  output {
    File? scores = "~{basename}.~{mode}.scores.hdf5"
    File? annots = "~{basename}.~{mode}.annot.hdf5"
    File output_vcf = "~{basename}.~{mode}.vcf.gz"
    File output_vcf_index = "~{basename}.~{mode}.vcf.gz.tbi"
  }
  runtime {
    docker: gatk_docker
    disks: "local-disk " + disk_size + " LOCAL"
    memory: memory_mb + " MiB"
  }
}

task SplitIntervalList {
  input {
    File interval_list
    Int scatter_count
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    Int disk_size
    String scatter_mode = "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.6.1"
  }

  parameter_meta {
    interval_list: {
      localization_optional: true
    }
  }

  command <<<
    gatk --java-options "-Xms3000m -Xmx3250m" SplitIntervals \
      -L ~{interval_list} -O  scatterDir -scatter ~{scatter_count} -R ~{ref_fasta} \
      -mode ~{scatter_mode} --interval-merging-rule OVERLAPPING_ONLY
    >>>

  runtime {
    memory: "3750 MiB"
    preemptible: 1
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " LOCAL"
    docker: gatk_docker
  }

  output {
    Array[File] output_intervals = glob("scatterDir/*")
  }
}

task SelectVariants {
  input {
    Array[File] input_intervals
    File input_vcf_file
    File input_vcf_index_file

    File ref_fasta
    File ref_fasta_index
    File ref_dict
  }

  String base_vcf = basename(input_vcf_file)
  Boolean is_compressed = basename(base_vcf, "gz") != base_vcf
  String vcf_suffix = if is_compressed then ".vcf.gz" else ".vcf"
  String vcf_index_suffix = if is_compressed then ".tbi" else ".idx"
  Int disk_size = 3 * ceil(size(input_vcf_file, "GiB") + size(input_vcf_index_file, "GiB") + size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB"))

  command <<<
    set -eo pipefail
    cat > /tmp/script <<END
%s
    END
    python3 /tmp/script
  >>>

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.2.4.1"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    memory: "3500 MiB"
    preemptible: 1
  }

  output {
    Array[File] output_vcf = glob("vcf_splits/*.${vcf_suffix}")
    Array[File] output_vcf_index = glob("vcf_splits/*.${vcf_index_suffix}")
  }
}