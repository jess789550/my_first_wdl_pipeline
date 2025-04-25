version 1.0

workflow gatk {
    input {
        File input_bam
        String file_prefix
        File ref_genome
        File ref_genome_bwt
        File ref_genome_sa
        File ref_genome_amb
        File ref_genome_ann
        File ref_genome_pac
    }

    call gatk_haplotypecaller {
        input:
            input_bam = input_bam,
            file_prefix = file_prefix,
            ref_genome = ref_genome,
            ref_genome_bwt = ref_genome_bwt,
            ref_genome_sa = ref_genome_sa,
            ref_genome_amb = ref_genome_amb,
            ref_genome_ann = ref_genome_ann,
            ref_genome_pac = ref_genome_pac
    }

    output {
        File vcf_file = gatk_haplotypecaller.vcf_file
    }
}

task gatk_haplotypecaller {
    input {
        File input_bam
        String file_prefix
        File ref_genome
        File ref_genome_bwt
        File ref_genome_sa
        File ref_genome_amb
        File ref_genome_ann
        File ref_genome_pac
    }
    
    command {
        /data/resources/GenomeAnalysisTK-4.5.0.0/gatk-4.5.0.0/gatk21 --java-options "-Xmx4g" HaplotypeCaller  \
        -R ${ref_genome} \
        -I ${input_bam} \
        -O ${file_prefix}.vcf.gz \
        -ERC GVCF 
    }
    
    output {
        File vcf_file = "${file_prefix}.vcf.gz"
    }
}
