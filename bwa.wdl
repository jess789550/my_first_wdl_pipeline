version 1.0

workflow bwa {
    input {
        File input_fastq_1
        File input_fastq_2
        String file_prefix
        File ref_genome
        File ref_genome_bwt
        File ref_genome_sa
        File ref_genome_amb
        File ref_genome_ann
        File ref_genome_pac
    }
    
    call bwa_mem {
        input:
            input_fastq_1 = input_fastq_1,
            input_fastq_2 = input_fastq_2,
            file_prefix = file_prefix,
            ref_genome = ref_genome,
            ref_genome_bwt = ref_genome_bwt,
            ref_genome_sa = ref_genome_sa,
            ref_genome_amb = ref_genome_amb,
            ref_genome_ann = ref_genome_ann,
            ref_genome_pac = ref_genome_pac
    }

    call sam_to_bam {
        input:
            sam_file = bwa_mem.file1,
            file_prefix = file_prefix
    }

    call sort_bam {
        input:
            bam_file = sam_to_bam.file2,
            file_prefix = file_prefix
    }

    call index_bam {
        input:
            sort_bam_file = sort_bam.sorted_bam,
            file_prefix = file_prefix
    }

    output {
        File sorted_bam = sort_bam.sorted_bam
        File sorted_bam_index = index_bam.bam_index
    }
}

task bwa_mem {
    input {
        File input_fastq_1
        File input_fastq_2
        String file_prefix
        File ref_genome
        File ref_genome_bwt
        File ref_genome_sa
        File ref_genome_amb
        File ref_genome_ann
        File ref_genome_pac
    }
    
    command {
        bwa mem -t 6 -R '@RG\tID:sample1\tSM:sample1\tPL:ILLUMINA\tLB:lib1\tPU:unit1' ${ref_genome} ${input_fastq_1} ${input_fastq_2} > ${file_prefix}.sam
    }
    
    output {
        File file1 = "${file_prefix}.sam"
    }
}

task sam_to_bam {
    input {
        File sam_file
        String file_prefix
    }
    
    command {
        samtools view -b -S ${sam_file} > ${file_prefix}.bam
    }
    
    output {
        File file2 = "${file_prefix}.bam"
    }
}

task sort_bam {
    input {
        File bam_file
        String file_prefix
    }
    
    command {
        samtools sort -o ${file_prefix}-sorted.bam ${bam_file}
    }
    
    output {
        File sorted_bam = "${file_prefix}-sorted.bam"
    }
}

task index_bam {
    input {
        File sort_bam_file
        String file_prefix
    }
    
    command {
        samtools index ${sort_bam_file} ${file_prefix}-sorted.bam.bai
    }
    
    output {
        File bam_index = "${file_prefix}-sorted.bam.bai"
    }
}
