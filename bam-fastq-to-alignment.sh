#!/usr/bin/env bash
set -euo pipefail

### ====== CONFIG ====== ###
BAM_DIR="/home/sana/Documents/GeneticCodon/collaborations/Medullablastoma_project/medullablastoma/sc-RNAseq/BAM"                          # folder with *.bam
BAMTOFASTQ="${BAM_DIR}/bamtofastq_linux"         # path to the binary you downloaded
FASTQ_DIR="${BAM_DIR}/../fastqs"                 # output FASTQs (one subfolder per sample)
STAR_INDEX="/home/sana/genome_reference/STAR_index"       # your GRCh38 STAR index
WHITELIST="/home/sana/Documents/whitelist/3M-february-2018.txt"    # 10x v3 whitelist
STAR_OUT="${BAM_DIR}/../star_out"                # STARsolo outputs
THREADS=12                                       # adjust to your CPU
BTQ_THREADS=8                                    # threads for bamtofastq (I/O heavy, not worth maxing)
### ===================== ###

mkdir -p "${FASTQ_DIR}" "${STAR_OUT}"

check_tools() {
  command -v STAR >/dev/null 2>&1 || { echo "STAR not found in PATH"; exit 1; }
  [[ -x "${BAMTOFASTQ}" ]] || { echo "bamtofastq binary missing or not executable: ${BAMTOFASTQ}"; exit 1; }
  [[ -d "${STAR_INDEX}" ]] || { echo "STAR index not found: ${STAR_INDEX}"; exit 1; }
  [[ -f "${WHITELIST}" ]] || { echo "Whitelist not found: ${WHITELIST}"; exit 1; }
}

convert_one_bam() {
  local bam="$1"
  local base
  base="$(basename "${bam}" .bam)"
  local fqdir="${FASTQ_DIR}/${base}"

  # If merged FASTQs already exist, skip conversion
  if [[ -f "${fqdir}/${base}_R1.fastq.gz" && -f "${fqdir}/${base}_R2.fastq.gz" ]]; then
    echo "[SKIP] FASTQs already present for ${base}"
    return 0
  fi

  # If output dir exists (from a previous partial run), move it aside
  if [[ -d "${fqdir}" ]]; then
    local bak="${fqdir}_old_$(date +%Y%m%d-%H%M%S)"
    echo "[INFO] ${fqdir} exists; moving to ${bak}"
    mv "${fqdir}" "${bak}"
  fi

  echo "[BTQ] Converting ${base} → ${fqdir}"
  # DO NOT pre-create ${fqdir}; bamtofastq wants to create it itself
  "${BAMTOFASTQ}" --nthreads="${BTQ_THREADS}" "${bam}" "${fqdir}/"

  echo "[MERGE] Merging chunks for ${base}"
  # Merge all parts per read type into single files
  if compgen -G "${fqdir}"/*_I1_*.fastq.gz > /dev/null; then
    cat "${fqdir}"/*_I1_*.fastq.gz > "${fqdir}/${base}_I1.fastq.gz"
  fi
  cat "${fqdir}"/*_R1_*.fastq.gz > "${fqdir}/${base}_R1.fastq.gz"
  cat "${fqdir}"/*_R2_*.fastq.gz > "${fqdir}/${base}_R2.fastq.gz"

  # Remove chunk files (comment out if you want to keep them)
  find "${fqdir}" -maxdepth 1 -type f -name "*_I1_*.fastq.gz" -delete 2>/dev/null || true
  find "${fqdir}" -maxdepth 1 -type f -name "*_R1_*.fastq.gz" -delete
  find "${fqdir}" -maxdepth 1 -type f -name "*_R2_*.fastq.gz" -delete

  # Quick sanity check: R1 and R2 line counts
  local r1 r2
  r1=$(zcat "${fqdir}/${base}_R1.fastq.gz" | wc -l)
  r2=$(zcat "${fqdir}/${base}_R2.fastq.gz" | wc -l)
  if [[ $((r1%4)) -ne 0 || $((r2%4)) -ne 0 || "${r1}" -ne "${r2}" ]]; then
    echo "[WARN] Read count mismatch for ${base}: R1=${r1} lines, R2=${r2} lines"
  else
    echo "[OK] ${base} FASTQs look consistent"
  fi
}

starsolo_one_sample() {
  local base="$1"
  local fqdir="${FASTQ_DIR}/${base}"
  local outprefix="${STAR_OUT}/${base}_"

  if [[ -d "${STAR_OUT}" && -f "${outprefix}Log.final.out" ]]; then
    echo "[SKIP] STARsolo results already present for ${base}"
    return 0
  fi

  echo "[STAR] Running STARsolo for ${base}"
  # Important: R2 (cDNA) first, R1 (CB+UMI) second
  STAR \
    --genomeDir "${STAR_INDEX}" \
    --readFilesIn "${fqdir}/${base}_R2.fastq.gz" "${fqdir}/${base}_R1.fastq.gz" \
    --readFilesCommand zcat \
    --runThreadN "${THREADS}" \
    --soloType Droplet \
    --soloCBwhitelist "${WHITELIST}" \
    --soloFeatures Gene \
    --soloUMIfiltering MultiGeneUMI \
    --soloCellFilter EmptyDrops_CR \
    --outFileNamePrefix "${outprefix}"
}

main() {
  check_tools

  shopt -s nullglob
  mapfile -t bams < <(ls "${BAM_DIR}"/*.bam)
  if [[ "${#bams[@]}" -eq 0 ]]; then
    echo "No BAM files found in ${BAM_DIR}"
    exit 1
  fi

  echo "[INFO] Found ${#bams[@]} BAMs"

  # 1) BAM -> FASTQ (+merge)
  for bam in "${bams[@]}"; do
    convert_one_bam "${bam}"
  done

  # 2) STARsolo for each sample
  for bam in "${bams[@]}"; do
    base="$(basename "${bam}" .bam)"
    starsolo_one_sample "${base}"
  done

  echo "[DONE] All samples processed."
}

main "$@"
