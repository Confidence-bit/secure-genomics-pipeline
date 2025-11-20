# Dockerfile - secure-genomics-pipeline (clean)
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Install system deps (may require > container privileges)
RUN apt-get update -y && apt-get install -y --no-install-recommends \
    python3 python3-pip python3-venv wget curl gzip nano ca-certificates \
    bwa samtools bcftools fastp mafft ivar \
    git build-essential \
  && apt-get clean && rm -rf /var/lib/apt/lists/*

# Python deps for the script
RUN pip3 install --no-cache-dir cryptography pyyaml

WORKDIR /pipeline

# Copy pipeline into container
COPY . /pipeline

# Ensure scripts are executable 
RUN chmod +x /pipeline/run_pipeline.py || true

# Default command: run pipeline 
CMD ["python3", "/pipeline/run_pipeline.py"]
