version 1.0

workflow vep {
    input {
        File vcf
        String file_prefix
    }

    call vep_annotation {
        input:
            vcf = vcf,
            file_prefix = file_prefix
    }

    output {
        File annotated_file = vep_annotation.annotated_file
    }
}

task vep_annotation {
    input {
        File vcf
        String file_prefix
    }
    
    command {
        /data/software_tmp/ensembl-vep-release-108/vep -i ${vcf} -o ${file_prefix}.txt --cache --dir_cache /data/vep/
    }
    
    output {
        File annotated_file = "${file_prefix}.txt"
    }
}
