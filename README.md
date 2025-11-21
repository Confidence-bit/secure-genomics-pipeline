# secure-genomics-pipeline
A beginner-friendly, reproducible, and secure bioinformatics workflow for downloading, processing, aligning, and protecting genomic sequence data.

**Project Overview**
This project demonstrates how to build a secure and reproducible genomic data processing pipeline using Python and basic command-line tools.
It highlights practical skills in:

**Bioinformatics Capabilities** 

This pipeline follows the major steps of a real-world genomics workflow:

**1. Sequence Acquisition**

  •	Downloads FASTA/FASTQ files from public databases.

  •	Saves them in the raw_reads/ directory

**2. Quality Control & Trimming**

  •	Cleans low-quality or invalid bases.

  •	Outputs saved in trimmed/.

**3. Reference Indexing**

  •	Prepares reference genome for alignment.

  •	Reference files stored in reference/.

**4. Sequence Alignment**
  •	Aligns sequences to the reference.

  •	SAM/BAM files generated in alignments

**5. SAM/BAM Processing and Indexing**

  •	Sorts, indexes, and processes alignment files. /.
  
**6. Multiple Sequence Alignment (MSA)**

  •	Compares sequences across samples.

  •	Stored in msa/.

**7. Consensus Sequence Generation**

  •	Builds consensus genomes from aligned reads.

  •	Saved in consensus/.

**8. Variant Calling**

  •	Detects SNPs, mutations, and genomic differences.

  •	VCF files written to variants/.

**Data Security (Two Encryption Methods)**

The pipeline includes two independent encryption approaches, demonstrating practical genomic data protection:

1. **OpenSSL Encryption (Command-line Level)**

  •	Used to encrypt sensitive variant files such as VCF.

  •	Produces encrypted files stored in encrypted/.

  •	Shows familiarity with standard cryptographic tools used in bioinformatics workflows.

2. **Python AES-256 Encryption (Fernet)**

  •	Encrypts selected output files using Python’s cryptography library.

  •	Demonstrates secure coding practices in Python.

  •	Keys handled separately to avoid plaintext exposure.

**Workflow Automation**

  •	All major steps are handled through Python scripts located in scripts/.

  •	Includes logging, validation, error handling, and step sequencing.

  •	The main controller script (run_pipeline.py) shows workflow orchestration.

**Reproducibility**

To ensure reproducible and transparent analysis:

  •	Structured directory layout for raw, processed, and encrypted data.

  •	Logging system (written to logs/) to track all steps, errors, and outputs.

  •	config/ folder contains parameters and settings used across the workflow.

Dockerfile is included for fully reproducible environments, useful for sharing or deploying the pipeline.

**Secure Coding Practices**

This project demonstrates:

  •	Input validation before processing data

  •	Avoiding plaintext storage of sensitive genomic results

  •	Separation of raw vs. processed vs. encrypted outputs

  •	Avoiding hard-coded keys

  •	Logging for traceability and auditability

  •	Correct use of cryptographic libraries and OpenSSL

  •	Git version control ensures all changes are tracked.

**Repository Structure**

secure-genomics-pipeline/

│── raw_reads/        # Downloaded FASTA/FASTQ files

│── trimmed/          # Cleaned/filtered sequences

│── reference/        # Reference genomes

│── alignments/       # Alignment outputs (SAM/FASTA, etc.)

│── consensus/        # Consensus sequences

│── msa/              # Multiple sequence alignment files

│── variants/         # Variant calling results (VCF)

│── encrypted/        # AES-encrypted output files

│── logs/             # Pipeline logs for auditing

│── scripts/          # Python workflow scripts

│── config/           # YAML or JSON configuration files

│── Dockerfile        # Containerized workflow definition

│── README.md         # Project documentation

**Features**

**1. Automated Genomic Data Download**

Fetches open-source sequences (e.g., SARS-CoV-2 FASTA files).

**2. Data Cleaning and Validation**

Ensures high-quality and safe-to-process nucleotide sequences.

**3. Alignment, MSA, Consensus, and Variant Calling**

Shows end-to-end genomic data analysis.

**4. Secure Logging**

Everything is logged for transparency, debugging, and audit trails.

**5. Dual Encryption System**

  •	OpenSSL for VCF security

  •	Python AES/Fernet for file-level encryption

**6. Dockerized Workflow**

Ensures reproducibility across systems.

**How This Workflow Ensures Security, Integrity & Reproducibility**

1. **Reproducible Workflow Design**

The pipeline follows a clear step-by-step process (download → clean → align → encrypt). Each step is implemented as a standalone script, ensuring that the same input always produces the same output, making the results easy to reproduce on any machine.

**Key features:**

  •	Consistent folder structure

  •	Clear input/output dependencies

  •	Pipeline logic contained in a single run_pipeline.py 

This ensures that any researcher—or security analyst—can reproduce the exact results on any machine.

2. **Version Control for Full Traceability**

All code, configurations, and documentation are managed using Git and GitHub.

This allows:

  •	Tracking of every change made to the workflow

  •	Rollback to previous versions

  •	Collaboration without losing project history

Version control provides a verifiable audit trail, which is essential for secure bioinformatics workflows.

3. **Secure Logging & Audit Trail**

The pipeline uses Python’s logging module to record:

  •	Which data was processed

  •	When each step ran

  •	Success/failure events

Logs create an audit trail, which is critical in secure genomics environments where data must be monitored, verified, and traceable.

4. **Data Encryption for Confidential Genomic Files**

Sensitive files (FASTA/FASTQ, alignment results, variants) can be encrypted using AES-256 through Python’s cryptography library.

This protects:

  •	Raw genomic reads

  •	Processed sequences

  •	Results files

  •	Keys

Encryption ensures genomic data remains protected at rest and cannot be accessed without the correct decryption key.

5. **Clean Separation of Data & Code**

Large or private genomic files are excluded from Git using a .gitignore.

This prevents accidental upload of:

  •	Sequencing data

  •	Variant files

  •	BAM/SAM/VCF

  •	Encryption keys

This follows real-world data governance principles, keeping the repository safe.

6. **Docker for Environment Reproducibility**
   
The project includes a Dockerfile (or supports Docker) that standardizes:

  •	Operating system

  •	Python version

  •	Dependencies

  •	Pipeline environment

This ensures the workflow runs identically on:

  •	Linux

  •	Windows

  •	macOS

  •	Cloud servers

Containerization eliminates “it works on my machine” problems and supports long-term reproducibility. Thereby ensuring data integrity, security, and reproducibility,

**How to Run the Pipeline**

**1. Clone the Repository**

git clone https://github.com/Confidence-bit/secure-genomics-pipeline.git

cd secure-genomics-pipeline

**2. Run the Pipeline Locally**

python3 run_pipeline.py

**3. (Optional) Run with Docker**

docker build -t secure-pipeline .

docker run secure-pipeline

**Security Features**

•	AES-256 encryption for sensitive data

•	Secure audit logs

•	Structured directory separation

•	Reproducible container image

**Author**

**Orji Confidence Ogechi**

GitHub: https://github.com/Confidence-bit






