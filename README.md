# my_first_wdl_pipeline
Building a simple pipeline using WDL

## Set up
- Download reference genome from: https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/
- Download FASTQ files from: https://github.com/hartwigmedical/testdata (I used 100K HiSeq)
- Download latest Cromwell jar files from: https://github.com/broadinstitute/cromwell/releases
- Download latest DNAnexus dxCompiler jar file from: https://github.com/dnanexus/dxCompiler/releases
- Create virtual environment with Java 11+, Docker, and DNAnexus installed

```
# Index reference genome file
bwa index Homo_sapiens.GRCh38.dna.alt.fa

samtools faindex Homo_sapiens.GRCh38.dna.alt.fa

java -jar /opt/software/picard-tools-1.96/CreateSequenceDictionary.jar R=Homo_sapiens.GRCh38.dna.alt.fa O=Homo_sapiens.GRCh38.dna.alt.dict
```

## Run pipeline locally
```
java -jar cromwell-88.jar run pipeline.wdl -i pipeline.json -o options.json  # Requires Java 11+
```

## Run pipeline on DNAnexus
```
dx login

dx upload Homo_Sapiens*

dx upload TESTX*

dx upload *.json

java -jar dxCompiler-2.13.0.jar compile pipeline.wdl -inputs inputs.json --project Jess-training-resources --destination /test_pipeline  # Requires Java 8 or 11

dx run pipeline -f inputs.dx.json
```
