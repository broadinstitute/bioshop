version 1.0

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

task bioshop_variant_fit_task {
  input {
    Array[File] classifier_training_data
    Boolean classifier_model_combine
    String classifier_model_type

    Int disk_size
    String bioshop_docker = "gcr.io/broad-jukebox-vm/bioshop:latest"
    String output_basename
  }

  String output_path = output_basename + "_" + classifier_model_type + "-classifier.pickle"
  String combine_arg = if classifier_model_combine then "--combine" else ""

  command <<<
    set -e -o pipefail

    newt fit \
      ~{combine_arg} \
      --classifier ~{classifier_model_type} \
      ~{sep=" " prefix("-i ", classifier_training_data)} \
      -o ~{output_path}
    >>>

  runtime {
    memory: "64G"
    cpu: "16"
    preemptible: 0
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    docker: bioshop_docker
  }

  output {
    File classifier_path = output_path
  }
}

task bioshop_variant_etl_task {
  input {
    File interval_list
    File query_vcf
    File query_vcf_index
    File target_vcf
    File target_vcf_index
    Array[String] interval_strats
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    Int disk_size
    String bioshop_docker = "gcr.io/broad-jukebox-vm/bioshop:latest"
    String output_basename
  }

  String etl_dataframe = output_basename + "-etl-df.pickle"

  command <<<
    set -e -o pipefail

    touch ~{query_vcf_index} ~{target_vcf_index}

    newt etl \
      --query_vcf ~{query_vcf} \
      --query_vcf_index ~{query_vcf_index} \
      --target_vcf ~{target_vcf} \
      --target_vcf_index ~{target_vcf_index} \
      --reference ~{ref_fasta} \
      --reference_index ~{ref_fasta_index} \
      -I ~{interval_list} \
      ~{sep=" " interval_strats} \
      -o ~{etl_dataframe}
    >>>

  runtime {
    memory: "48G"
    cpu: "16"
    preemptible: 0
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " LOCAL"
    docker: bioshop_docker
  }

  output {
    File etl_dataframe = etl_dataframe
  }
}

workflow bioshop_variant_etl {
    String pipeline_version = "1.0.0"

    input {
        File unpadded_intervals
        Int scatter_count = 10

        File query_vcf
        File query_vcf_index
        File target_vcf
        File target_vcf_index

        File ref_fasta
        File ref_fasta_index
        File ref_dict
        Array[Pair[String, String]] interval_strat_map

        Boolean classifier_model_combine = true
        String classifier_model_type = "xgb"

        Int small_disk = 128
        Int medium_disk = 256
        Int large_disk = 1024
        Int huge_disk = 2048
    }

    String output_basename = basename(query_vcf)

    call SplitIntervalList as interval_shards {
      input:
        interval_list = unpadded_intervals,
        scatter_count = scatter_count,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        disk_size = small_disk,
    }

    scatter (pair in interval_strat_map) {
      String interval_strats = "-S ~{pair.left}=~{pair.right}"
    }
    
    scatter (interval_list in interval_shards.output_intervals) {
      call bioshop_variant_etl_task as etl {
        input:
          interval_list = interval_list,
          query_vcf = query_vcf,
          query_vcf_index = query_vcf_index,
          target_vcf = target_vcf,
          target_vcf_index = target_vcf_index,
          interval_strats = interval_strats,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          disk_size = large_disk,
          output_basename = output_basename
      }
    }

    call bioshop_variant_fit_task as fit {
      input:
        classifier_training_data = etl.etl_dataframe,
        classifier_model_type = classifier_model_type,
        classifier_model_combine = classifier_model_combine,
        disk_size = large_disk,
        output_basename = output_basename
    }

    output {
        File classifier_path = fit.classifier_path
    }
}
