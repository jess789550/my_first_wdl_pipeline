version 1.0

workflow gatk {
    input {
        File input_bam
        File bai_file
        String file_prefix
        File ref_genome
        File ref_genome_fai
        File ref_genome_dict
    }

    call gatk_haplotypecaller {
        input:
            input_bam = input_bam,
            bai_file = bai_file,
            file_prefix = file_prefix,
            ref_genome = ref_genome,
            ref_genome_fai = ref_genome_fai,
            ref_genome_dict = ref_genome_dict
    }

    output {
        File vcf_file = gatk_haplotypecaller.vcf_file
    }
}

task gatk_haplotypecaller {
    input {
        File input_bam
        File bai_file
        String file_prefix
        File ref_genome
        File ref_genome_fai
        File ref_genome_dict
    }
    
    command {
        cp ${input_bam} ./input.bam
        cp ${bai_file} ./input.bam.bai
        
        gatk --java-options "-Xmx4g" HaplotypeCaller  \
        -R ${ref_genome} \
        -I input.bam \
        -O ${file_prefix}.vcf.gz 
    }
    
    output {
        File vcf_file = "${file_prefix}.vcf.gz"
    }
    
    runtime {
        docker: "broadinstitute/gatk:4.5.0.0"
    }
}
