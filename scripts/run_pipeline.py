import subprocess
import yaml
from logging_module import log
from encrypt import encrypt_file

# Load config.yaml
with open("config/config.yaml") as f:
    config = yaml.safe_load(f)

REF = config["reference"]
SAMPLE = config["sample_name"]

RAW_R1 = "raw_reads/sample_R1.fastq.gz"
RAW_R2 = "raw_reads/sample_R2.fastq.gz"

TRIM_R1 = "trimmed/sample_R1_trimmed.fastq.gz"
TRIM_R2 = "trimmed/sample_R2_trimmed.fastq.gz"

BAM = "alignments/sample.sorted.bam"
CONS = "consensus/sample.fa"
MSA = "msa/aligned.fasta"
VCF = "variants/sample.vcf.gz"
ENC = "encrypted/sample.enc"

def run(cmd, step):
    log(f"START: {step}")
    subprocess.run(cmd, shell=True, check=True)
    log(f"DONE: {step}")

# Pipeline
run(f"fastp -i {RAW_R1} -I {RAW_R2} -o {TRIM_R1} -O {TRIM_R2}", "QC + TRIMMING")
run(f"bwa mem {REF} {TRIM_R1} {TRIM_R2} | samtools sort -o {BAM}", "ALIGNMENT")
run(f"samtools index {BAM}", "BAM INDEX")
run(f"samtools mpileup -A -d 0 -Q 0 -f {REF} {BAM} | ivar consensus -p consensus/sample", "CONSENSUS")
run(f"mafft --auto reference/reference.fasta consensus/sample.fa > {MSA}", "MSA")
run(f"bcftools mpileup -f {REF} {BAM} | bcftools call -mv -Oz -o {VCF}", "VARIANT CALLING")

encrypt_file(VCF, ENC)
log("Pipeline complete.")
