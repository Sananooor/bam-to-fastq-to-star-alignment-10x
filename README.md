# BAM to FASTQ to STARsolo Alignment (10x Genomics scRNA-seq)

This repository provides a robust, resume-safe, and parallelised Bash script for converting BAM files (from 10x Genomics single-cell RNA-seq) to FASTQ, and then aligning them with [STARsolo](https://github.com/alexdobin/STAR). The workflow handles multiple samples in parallel, merges output chunks, cleans up intermediary files, and is safe to rerun.

## Features

- **Parallelised BAM→FASTQ conversion and STARsolo alignment** with job limits to avoid overloading your system.
- **Resume-safe:** Skips samples that are already processed, allowing you to rerun the script without redoing completed work.
- **Chunk merging:** Automatically merges per-sample FASTQ chunks and cleans up chunk files.
- **Sanity checks:** Verifies matching read pairs and warns about inconsistencies.
- **Customizable paths and resource limits.**

## Workflow Overview

1. **BAM to FASTQ:** Converts each BAM to 10x-style FASTQ files using `bamtofastq` (Cell Ranger binary or similar).
2. **FASTQ Merging:** Merges chunked FASTQ files into final R1, R2, and I1 files per sample.
3. **Alignment:** Runs STARsolo for each sample, producing gene count matrices and reports.

---

## Requirements

- **[STAR](https://github.com/alexdobin/STAR)** (in your `$PATH`)
- **bamtofastq** binary (downloaded from 10x Genomics or Cell Ranger)
- **BAM files** produced by 10x Genomics pipelines
- **STAR genome index** for your reference genome
- **Cell barcode whitelist** (e.g., `3M-february-2018.txt` for 10x v3 chemistry)
- **Bash** (script tested with bash 4+)
- Sufficient storage and compute resources

---

## Quick Start

1. **Clone this repository** and make the script executable:
   ```bash
   git clone https://github.com/Sananooor/bam-to-fastq-to-star-alignment-10x.git
   cd bam-to-fastq-to-star-alignment-10x
   chmod +x bam-fastq-to-alignment.sh
   ```

2. **Edit Configuration**  
   Open `bam-fastq-to-alignment.sh` and update the following variables in the CONFIG section:
   ```bash
   BAM_DIR="~/sc-RNAseq/BAM"                 # Directory with input BAM files
   BAMTOFASTQ="${BAM_DIR}/bamtofastq_linux"  # Path to bamtofastq binary
   FASTQ_DIR="${BAM_DIR}/../fastqs"          # Output dir for FASTQs
   STAR_INDEX="/path/to/STAR_index"          # STAR genome index
   WHITELIST="/path/to/whitelist.txt"        # Barcode whitelist (10x)
   STAR_OUT="${BAM_DIR}/../star_out"         # STARsolo output dir

   THREADS=12          # Threads for STAR
   BTQ_THREADS=8       # Threads for bamtofastq
   CONVERT_JOBS=2      # BAM→FASTQ jobs in parallel
   ALIGN_JOBS=2        # STARsolo jobs in parallel
   ```

3. **Prepare Directories**  
   Make sure all paths exist and the required files are accessible.

4. **Run the Pipeline**
   ```bash
   ./bam-fastq-to-alignment.sh
   ```

---

## Output

- **FASTQ files:** Saved per-sample in `${FASTQ_DIR}/SAMPLE_NAME/`
- **STARsolo results:** Per-sample output in `${STAR_OUT}/`
- **Logs:** Progress and warnings are printed to the terminal

---

## Tips and Notes

- **Resume-safety:** You can safely rerun the script; it will skip completed samples.
- **Partial runs:** If a sample partially failed, its directory is backed up and reprocessed.
- **Resource tuning:** Adjust `THREADS`, `BTQ_THREADS`, `CONVERT_JOBS`, and `ALIGN_JOBS` to match your hardware.
- **10x Chemistry:** Make sure your whitelist matches your chemistry version (v2/v3).

---

## Troubleshooting

- **Missing tools:** The script will check for `STAR`, `bamtofastq`, STAR index, and whitelist before running.
- **Performance:** For large datasets, ensure you have enough CPU, memory, and disk space.
- **Output files:** If output files are not generated, check for errors in the terminal and ensure all paths are correct.

---

## License

This script is provided as-is for research purposes. Please cite the original tools (**STAR**, **bamtofastq/Cell Ranger**) as appropriate.

---

## Author

[Sananooor](https://github.com/Sananooor)

---

**Contributions and issues are welcome!**
