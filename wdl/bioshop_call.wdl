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

task gather_vcfs {
  input {
    Array[File] input_vcfs
    String output_vcf_basename

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.6.1"
    Int cpu = 1
    Int memory_mb = 16000
    Int disk_size = disk_size
  }
  Int command_mem = memory_mb - 1000
  Int max_heap = memory_mb - 500

  command <<<
    set -e -o pipefail

    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
    GatherVcfs \
    -I ~{sep=' -I ' input_vcfs} \
    --REORDER_INPUT_BY_FIRST_VARIANT \
    -O ~{output_vcf_basename}.vcf.gz

    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
    IndexFeatureFile -I ~{output_vcf_basename}.vcf.gz
  >>>
  runtime {
    docker: gatk_docker
    disks: "local-disk ${disk_size} LOCAL"
    memory: "${memory_mb} MiB"
    cpu: cpu
  }
  output {
    File output_vcf = "~{output_vcf_basename}.vcf.gz"
    File output_vcf_index = "~{output_vcf_basename}.vcf.gz.tbi"
  }
}

task bioshop_variant_call_task {
  input {
    File interval_list
    File query_vcf
    File query_vcf_index
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File model
    Array[String] interval_strats
    Int disk_size
    String bioshop_docker = "gcr.io/broad-jukebox-vm/bioshop:latest"
    String output_basename
  }

  String output_vcf = output_basename + "-called.vcf"

  command <<<
    set -e -o pipefail

    touch ~{query_vcf} ~{target_vcf}

    newt call \
      --query_vcf ~{query_vcf} \
      --query_vcf_index ~{query_vcf_index} \
      --reference ~{ref_fasta} \
      --reference_index ~{ref_fasta_index} \
      --classifier ~{model} \
      -I ~{interval_list} \
      ~{sep=" " interval_strats} \
      -o ~{output_vcf}
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
    File output_vcf = output_vcf
  }
}

workflow bioshop_variant_call {
    String pipeline_version = "1.0.0"

    input {
        File unpadded_intervals
        Int scatter_count = 10

        File query_vcf
        File query_vcf_index
        String model

        File ref_fasta
        File ref_fasta_index
        File ref_dict
        Array[Pair[String, String]] interval_strat_map


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
      call bioshop_variant_call_task as call_task {
        input:
          query_vcf = query_vcf,
          query_vcf_index = query_vcf_index,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          model = model,
          interval_strats = interval_strats,
          interval_list = interval_list,
          disk_size = large_disk,
          output_basename = output_basename
      }
    }

    call gather_vcfs as gather {
      input:
        input_vcfs = call_task.output_vcf,
        output_vcf_basename = output_basename,
        disk_size = huge_disk
    }

    output {
        File output_vcf = gather.output_vcf
        File output_vcf_index = gather.output_vcf_index
    }
}
