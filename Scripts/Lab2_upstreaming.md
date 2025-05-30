``` bash 

conda activate RNA-seq # write your own environent name

# To download the reference transcriptome file on Ensembl, we selected the cDNA fasta file 
# we fetched the file using wget 
wget https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# then index the fasta file
kallisto index -i HS_Kallisto_index RNA-seq_project/Ref/Homo_sapiens.GRCh38.cdna.all.fa.gz


# # # # # # # # # # # Quality Control of Fastq Files (Raw Data) # # # # # # # # # # # # 
# 1. Generating Fastqc Reprots with Fastqc 
fastqc *.fastq.gz -o fastqc_dir -t 8 # (*.fastq.gz) this will collectively input all the file that extension is .fastq.gz

# to collect the fastqc reprots into one agregated report, we used multiqc 
multiqc -d fastqc_dir 

# After you identify the problem of your data, you decide to proceed with which tool
# 3. Filtering Raw Quality Reads using Fastp
conda install fastp -c bioconda 

fastp -i input_file_forward.fastq.gz -o output_filtered_forward.fastq.gz -I input_file_reverse.fastq.gz -O output_filtered_reverse.fastq.gz -q 30 # for minimum Q score of 30

# Now, After cleaning data you are ready for mapping reads using Kallisto 
# 4. Mapping reads 
kallisto quant -i HS_Kallisto_index -o This_should_be_Sample_ID_dir input_file_reverse.fastq.gz output_filtered_reverse.fastq.gz
# Usage: kallisto quant [arguments] FASTQ-files
# Required arguments:
# -i, --index=STRING            Filename for the kallisto index to be used for
#                              quantification
# -o, --output-dir=STRING       Directory to write output to
```
