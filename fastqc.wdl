version 1.0

workflow fastqc {
    input {
        File input_fastq_1
        File input_fastq_2
    }
    
    call fastqc_task {
        input:
            input_fastq_1 = input_fastq_1,
            input_fastq_2 = input_fastq_2
    }

    output {
        File file1 = fastqc_task.file1
        File file2 = fastqc_task.file2
    }
}


task fastqc_task {
    input {
        File input_fastq_1
        File input_fastq_2
    }
    
    command {
        fastqc ${input_fastq_1} > ${basename(input_fastq_1, ".fastq.gz")}.html
        fastqc ${input_fastq_2} > ${basename(input_fastq_2, ".fastq.gz")}.html
    }
    
    output {
        File file1 = "${basename(input_fastq_1, ".fastq.gz")}.html"
        File file2 = "${basename(input_fastq_2, ".fastq.gz")}.html"
    }
    
    runtime {
        docker: "biocontainers/fastqc:v0.11.9_cv8"
    }
}

