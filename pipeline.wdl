# Run Cromwell workflow: java -jar cromwell-88.jar run pipeline.wdl -i pipeline.json -o options.json
# java -jar dxCompiler-2.13.0.jar compile pipeline.wdl -inputs inputs.json 
# ref_genome from https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/
# FASTQ files from https://github.com/hartwigmedical/testdata

version 1.0

# Import external WDL files
import "fastqc.wdl" as fastqc
import "bwa.wdl" as bwa
import "gatk.wdl" as gatk

workflow pipeline {
  input {
    # Input file for the reference genome and indices
    File ref_genome
    File ref_genome_bwt
    File ref_genome_sa
    File ref_genome_amb
    File ref_genome_ann
    File ref_genome_pac
    
    # Define input_fastq as an array of FilePair records
    Array[Pair[File,File]] input_fastq
  }
  
  # Scatter over pairs of input fastq files
  scatter (pair in input_fastq) {
    String file_prefix = basename(pair.left, ".fastq.gz")
    
    # Initialise task for each pair of fastq files
    call initialise {
      input: 
        input_fastq_1 = pair.left,
        input_fastq_2 = pair.right,
        ref_genome = ref_genome
    }
    
    # FastQC task for each pair of fastq files
    call fastqc.fastqc {
      input:
        input_fastq_1 = pair.left,
        input_fastq_2 = pair.right
    }
    
    # BWA MEM alignment
    call bwa.bwa {
      input:
        input_fastq_1 = pair.left,
        input_fastq_2 = pair.right,
        ref_genome = ref_genome,
        file_prefix = file_prefix,
        ref_genome_bwt = ref_genome_bwt,
        ref_genome_sa = ref_genome_sa,
        ref_genome_amb = ref_genome_amb,
        ref_genome_ann = ref_genome_ann,
        ref_genome_pac = ref_genome_pac
    }

    # GATK HaplotypeCaller Variant Calling
    call gatk.gatk_haplotypecaller {
      input:
        input_bam = bwa.file3,
        ref_genome = ref_genome,
        file_prefix = file_prefix,
        ref_genome_bwt = ref_genome_bwt,
        ref_genome_sa = ref_genome_sa,
        ref_genome_amb = ref_genome_amb,
        ref_genome_ann = ref_genome_ann,
        ref_genome_pac = ref_genome_pac
    }
    
    # Let user know pipeline is Finished
    call finish {
      input:
        gatk_output = gatk_haplotypecaller.vcf_file
    }
  }
}

task initialise {
  input {
    File input_fastq_1
    File input_fastq_2
    File ref_genome
  }

  command {
    echo 
      """ 
      -----------------
      Starting Pipeline
      -----------------   
      
      ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⢀⣠⣶⠾⠛⠛⠋⠉⠛⠛⠷⢶⣤⡀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⣠⣾⠟⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣷⣄⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⣴⠟⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣧⠀⠀⠀⠀⠀
⠀⠀⠀⠀⣸⡏⠀⠀⠀⠀⣴⣶⡆⠀⠀⠀⠀⢠⣶⣦⠀⠀⠀⠀⢹⣇⠀⠀⠀⠀
⠀⠀⠀⢠⣿⠀⠀⠀⠀⠀⠛⠟⠃⠀⠀⠀⠀⠘⠻⠛⠀⠀⠀⠀⠀⣿⡄⠀⠀⠀
⠀⠀⠀⢸⣿⠀⠀⠀⠀⢠⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡄⠀⠀⠀⠀⢸⡇⠀⠀⠀
⠀⠀⠀⠘⣿⠀⠀⠀⠀⠸⣧⠀⠀⠀⠀⠀⠀⠀⠀⣼⡏⠀⠀⠀⠀⣿⠇⠀⠀⠀
⠀⠀⠀⠀⢻⣇⠀⠀⠀⠀⠙⢷⣦⣀⡀⠀⣀⣠⡾⠋⠀⠀⠀⠀⢰⡟⠀⠀⠀⠀
⠀⠀⠀⠀⠀⢻⣦⠀⠀⠀⠀⠀⠈⠉⠙⠋⠉⠉⠀⠀⠀⠀⠀⣰⡟⠁⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠙⢷⣄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣠⣾⠏⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠙⠻⢶⣤⣄⣀⣀⣀⣀⣠⣤⡶⠟⠋⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
      """
    echo ${input_fastq_1}
    echo ${input_fastq_2}
    echo ${ref_genome}
  }

  output {
    String inputs = stdout()
  }
}

task finish {
  input{
    File gatk_output
  }

  command {
    echo
      """
      -----------------
      Finished pipeline
      -----------------

  ⠀⠀⠀⠀⠀⠀⠀⠀⢀⣠⣶⠾⠛⠛⠋⠉⠛⠛⠷⢶⣤⡀⠀⠀⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⠀⣠⣾⠟⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣷⣄⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⣴⠟⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣧⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⣸⡏⠀⠀⠀⠀⣴⣶⡆⠀⠀⠀⠀⢠⣶⣦⠀⠀⠀⠀⢹⣇⠀⠀⠀⠀
  ⠀⠀⠀⢠⣿⠀⠀⠀⠀⠀⠛⠟⠃⠀⠀⠀⠀⠘⠻⠛⠀⠀⠀⠀⠀⣿⡄⠀⠀⠀
  ⠀⠀⠀⢸⣿⠀⠀⠀⠀⢠⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡄⠀⠀⠀⠀⢸⡇⠀⠀⠀
  ⠀⠀⠀⠘⣿⠀⠀⠀⠀⠸⣧⠀⠀⠀⠀⠀⠀⠀⠀⣼⡏⠀⠀⠀⠀⣿⠇⠀⠀⠀
  ⠀⠀⠀⠀⢻⣇⠀⠀⠀⠀⠙⢷⣦⣀⡀⠀⣀⣠⡾⠋⠀⠀⠀⠀⢰⡟⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⢻⣦⠀⠀⠀⠀⠀⠈⠉⠙⠋⠉⠉⠀⠀⠀⠀⠀⣰⡟⠁⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⠀⠙⢷⣄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣠⣾⠏⠀⠀⠀⠀⠀⠀
  ⠀⠀⠀⠀⠀⠀⠀⠀⠙⠻⢶⣤⣄⣀⣀⣀⣀⣠⣤⡶⠟⠋⠀⠀⠀⠀⠀⠀⠀⠀
       
  ⠀⠀⠀⠀"""
  }
}
