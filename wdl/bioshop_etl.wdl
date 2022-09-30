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
    disks: "local-disk " + disk_size + " HDD"
    docker: gatk_docker
  }

  output {
    Array[File] output_intervals = glob("scatterDir/*")
  }
}

task BioshopETL {
  input {
    File interval_list
    File query_vcf
    File target_vcf
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
    newt etl \
      --query_vcf ~{query_vcf} \
      --target_vcf ~{target_vcf} \
      -I ~{interval_list} \
      ~{sep=" " interval_strats} \
      -o ~{etl_dataframe}
    >>>

  runtime {
    memory: "3750 MiB"
    preemptible: 1
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    docker: bioshop_docker
  }

  output {
    File etl_dataframe = etl_dataframe
  }
}

workflow VariantETL {
    String pipeline_version = "1.0.0"

    input {
        File unpadded_intervals
        Int scatter_count = 20

        File query_vcf
        File query_vcf_idx
        File target_vcf
        File target_vcf_idx

        File ref_fasta
        File ref_fasta_index
        File ref_dict
        Array[Pair[String, String]] interval_strat_map

        Int small_disk = 100
        Int medium_disk = 200
        Int large_disk = 1000
        Int huge_disk = 2000
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

    scatter (pair in interval_strat_map) {
      String interval_strats = "-S ~{pair.left}=~{pair.right}"
    }
    
    scatter (interval_list in interval_shards.output_intervals) {
      call BioshopETL {
        input:
          interval_list = interval_list,
          query_vcf = query_vcf,
          target_vcf = target_vcf,
          interval_strats = interval_strats,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          disk_size = 1000,
          output_basename = basename(query_vcf),
      }
    }

    output {
        Array[File] etl_dataframe = BioshopETL.etl_dataframe
    }
}
