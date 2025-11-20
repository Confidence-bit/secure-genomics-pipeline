#!/usr/bin/env python3
"""
run_pipeline.py
Robust single-file pipeline:
FASTQ -> fastp -> bwa -> samtools -> ivar consensus -> mafft MSA -> bcftools VCF -> encrypt (Fernet)
Requirements:
    python3, pyyaml, cryptography
Unix tools: fastp, bwa, samtools, ivar, mafft, bcftools
"""

import os
import sys
import subprocess
import datetime
from pathlib import Path

try:
    import yaml
except Exception:
    print("Missing Python dependency: pyyaml. Install with: pip3 install pyyaml", file=sys.stderr)
    sys.exit(1)

try:
    from cryptography.fernet import Fernet
except Exception:
    print("Missing Python dependency: cryptography. Install with: pip3 install cryptography", file=sys.stderr)
    sys.exit(1)


# ---------------------------
# Config (simple defaults; override in config/config.yaml)
# ---------------------------
CFG_PATH = Path("config/config.yaml")
DEFAULT = {
    "reference": "reference/reference.fasta",
    "sample_name": "sample",
    "log_file": "logs/pipeline.log"
}

if CFG_PATH.exists():
    with open(CFG_PATH) as f:
        cfg = yaml.safe_load(f) or {}
else:
    cfg = {}

REF = Path(cfg.get("reference", DEFAULT["reference"]))
SAMPLE = cfg.get("sample_name", DEFAULT["sample_name"])
LOGFILE = Path(cfg.get("log_file", DEFAULT["log_file"]))

# Paths (consistent names used in your runs)
RAW_R1 = Path(f"raw_reads/{SAMPLE}_R1.fastq.gz")
RAW_R2 = Path(f"raw_reads/{SAMPLE}_R2.fastq.gz")
TRIM_R1 = Path(f"trimmed/{SAMPLE}_R1_trimmed.fastq.gz")
TRIM_R2 = Path(f"trimmed/{SAMPLE}_R2_trimmed.fastq.gz")
SAM = Path(f"alignments/{SAMPLE}.sam")
BAM = Path(f"alignments/{SAMPLE}.sorted.bam")
CONSENSUS = Path(f"consensus/{SAMPLE}.fa")
MSA_DIR = Path("msa")
COMBINED = MSA_DIR / "combined.fasta"
ALIGNED = MSA_DIR / "aligned.fasta"
VCF = Path(f"variants/{SAMPLE}.vcf.gz")
ENCRYPTED = Path(f"encrypted/{SAMPLE}.vcf.gz.enc")
FERNET_KEY = Path("fernet.key")

# Ensure directories exist
for d in ["raw_reads", "trimmed", "reference", "alignments", "consensus", "msa", "variants", "encrypted", "logs"]:
    Path(d).mkdir(parents=True, exist_ok=True)


# ---------------------------
# Simple logging
# ---------------------------
def log(msg):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    try:
        LOGFILE.parent.mkdir(parents=True, exist_ok=True)
        with open(LOGFILE, "a") as fh:
            fh.write(line + "\n")
    except Exception as e:
        print(f"Warning: failed to write log: {e}", file=sys.stderr)


# ---------------------------
# Helpers
# ---------------------------
def run(cmd, step=None, check=True, env=None):
    if step:
        log(f"START: {step}")
    log(f"CMD: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=check, env=env)
    except subprocess.CalledProcessError as e:
        log(f"ERROR: command failed (exit {e.returncode}): {e}")
        raise
    if step:
        log(f"DONE: {step}")


def file_exists(path: Path, desc="file"):
    if not path.exists():
        raise FileNotFoundError(f"Required {desc} not found: {path}")


# ---------------------------
# Key functions
# ---------------------------
def ensure_reference_index(ref: Path):
    # samtools .fai
    if not ref.exists():
        raise FileNotFoundError(f"Reference fasta not found: {ref}")
    faidx = ref.with_suffix(ref.suffix + ".fai")
    if not faidx.exists():
        log("Creating samtools faidx for reference")
        run(f"samtools faidx {ref}", "samtools faidx")
    # bwa index: check for .bwt file
    bwt = ref.with_suffix(ref.suffix + ".bwt")
    if not bwt.exists():
        log("Creating BWA index for reference")
        run(f"bwa index {ref}", "bwa index")


def generate_fernet_key(keypath: Path):
    if keypath.exists():
        log(f"Fernet key already exists: {keypath}")
        return
    keypath.write_bytes(Fernet.generate_key())
    keypath.chmod(0o600)
    log(f"Fernet key generated at {keypath}")


def encrypt_with_fernet(keypath: Path, infile: Path, outfile: Path):
    file_exists(keypath, "fernet key")
    file_exists(infile, "input file to encrypt")
    key = keypath.read_bytes()
    f = Fernet(key)
    data = infile.read_bytes()
    out = f.encrypt(data)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    outfile.write_bytes(out)
    log(f"Encrypted {infile} -> {outfile}")


# ---------------------------
# Pipeline Steps (robust)
# ---------------------------
def main():
    try:
        # 0. sanity checks
        log("Pipeline starting")
        log(f"Sample: {SAMPLE}")
        log(f"Reference: {REF}")

        # 1. FASTQ files (raw)
        # Accept either sample_R1/sample_R2 (preferred) OR old names if present.
        if not RAW_R1.exists() or not RAW_R2.exists():
            # try alternate names (you used sample_R1.fastq.gz etc earlier)
            alt1 = Path(f"raw_reads/{SAMPLE}.fastq.gz")
            if alt1.exists():
                log(f"Found alternate raw file {alt1}; renaming to expected paired names is required.")
                raise FileNotFoundError("Paired raw FASTQ files not found. Ensure *_R1.fastq.gz and *_R2.fastq.gz exist in raw_reads/")
            else:
                raise FileNotFoundError(f"Paired raw FASTQ files not found: {RAW_R1}, {RAW_R2}")

        # 2. QC + trimming (fastp)
        trimmed_cmd = f"fastp -i {RAW_R1} -I {RAW_R2} -o {TRIM_R1} -O {TRIM_R2} --html trimmed/fastp_report.html --json trimmed/fastp_report.json"
        run(trimmed_cmd, "fastp QC + trimming")

        # 3. Reference checks & indexing
        ensure_reference_index(REF)

        # 4. Alignment (bwa mem) -> sam -> sorted bam
        run(f"bwa mem {REF} {TRIM_R1} {TRIM_R2} > {SAM}", "bwa mem")
        run(f"samtools view -bS {SAM} > alignments/{SAMPLE}.bam", "samtools view BAM")
        run(f"samtools sort alignments/{SAMPLE}.bam -o {BAM}", "samtools sort")
        run(f"samtools index {BAM}", "samtools index")

        # 5. Consensus (ivar)
        # Ensure ivar exists; iVar expects mpileup piped
        run(f"samtools mpileup -aa -A -d 0 -Q 0 -f {REF} {BAM} | ivar consensus -p consensus/{SAMPLE} -q 20 -t 0.6 -m 10", "iVar consensus")
        file_exists(CONSENSUS, "consensus FASTA")

        # 6. Prepare combined fasta for MSA (reference + consensus)
        MSA_DIR.mkdir(parents=True, exist_ok=True)
        log(f"Creating combined FASTA: {COMBINED}")
        with open(COMBINED, "w") as out:
            with open(REF, "r") as r:
                out.write(r.read().rstrip() + "\n")
            with open(CONSENSUS, "r") as c:
                out.write(c.read().rstrip() + "\n")
        log(f"Combined FASTA written to {COMBINED}")

        # 7. MSA (MAFFT)
        run(f"mafft --auto {COMBINED} > {ALIGNED}", "MAFFT alignment")
        file_exists(ALIGNED, "MAFFT output")

        # 8. Variant calling (bcftools)
        run(f"bcftools mpileup -f {REF} {BAM} | bcftools call -mv -Oz -o {VCF}", "bcftools call")
        run(f"bcftools index {VCF}", "bcftools index")
        log(f"VCF written: {VCF}")

        # 9. Encryption (Fernet Python)
        generate_fernet_key(FERNET_KEY)
        encrypt_with_fernet(FERNET_KEY, VCF, ENCRYPTED)

        log("Pipeline finished successfully")

    except Exception as e:
        log(f"Pipeline failed: {e}")
        raise


if __name__ == "__main__":
    main()
