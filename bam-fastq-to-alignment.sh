#!/usr/bin/env bash
# BAM → FASTQ → STARsolo (resume-safe + limited parallelism)
set -euo pipefail

### ========= CONFIG (EDIT ME) ========= ###
BAM_DIR="~/sc-RNAseq/BAM" # change the path accordingly
BAMTOFASTQ="${BAM_DIR}/bamtofastq_linux"        # downloaded binary
FASTQ_DIR="${BAM_DIR}/../fastqs"                # per-sample FASTQs
STAR_INDEX="/home/sana/genome_reference/STAR_index"
WHITELIST="/home/sana/Documents/whitelist/3M-february-2018.txt"   # 10x v3
STAR_OUT="${BAM_DIR}/../star_out"

THREADS=12          # threads for STAR
BTQ_THREADS=8       # threads for bamtofastq
CONVERT_JOBS=2      # parallel BAM→FASTQ jobs
ALIGN_JOBS=2        # parallel STARsolo jobs
### ==================================== ###

mkdir -p "${FASTQ_DIR}" "${STAR_OUT}"

check_tools() {
  command -v STAR >/dev/null 2>&1 || { echo "STAR not found in PATH"; exit 1; }
  [[ -x "${BAMTOFASTQ}" ]] || { echo "bamtofastq binary missing/not executable: ${BAMTOFASTQ}"; exit 1; }
  [[ -d "${STAR_INDEX}" ]] || { echo "STAR index not found: ${STAR_INDEX}"; exit 1; }
  [[ -f "${WHITELIST}" ]] || { echo "Whitelist not found: ${WHITELIST}"; exit 1; }
}

merge_fastqs() {  # merge_fastqs "<glob>" "<out>"
  local pattern="$1" out="$2"
  mapfile -t parts < <(compgen -G "$pattern") || true
  if ((${#parts[@]})); then
    cat "${parts[@]}" > "$out"
  else
    echo "[INFO] No matches for $pattern; skipping $out"
  fi
}

convert_one_bam() {
  local bam="$1"
  local base; base="$(basename "${bam}" .bam)"
  local fqdir="${FASTQ_DIR}/${base}"

  # already merged? skip
  if [[ -f "${fqdir}/${base}_R1.fastq.gz" && -f "${fqdir}/${base}_R2.fastq.gz" ]]; then
    echo "[SKIP] FASTQs already present for ${base}"
    return 0
  fi

  # if dir exists from a previous partial run, move aside
  if [[ -d "${fqdir}" ]]; then
    local bak="${fqdir}_old_$(date +%Y%m%d-%H%M%S)"
    echo "[INFO] ${fqdir} exists; moving to ${bak}"
    mv "${fqdir}" "${bak}"
  fi

  echo "[BTQ] Converting ${base} → ${fqdir}"
  # do NOT pre-create fqdir; bamtofastq needs to create it
  "${BAMTOFASTQ}" --nthreads="${BTQ_THREADS}" "${bam}" "${fqdir}/"

  echo "[MERGE] ${base}"
  merge_fastqs "${fqdir}/*_I1_*.fastq.gz" "${fqdir}/${base}_I1.fastq.gz"
  merge_fastqs "${fqdir}/*_R1_*.fastq.gz" "${fqdir}/${base}_R1.fastq.gz"
  merge_fastqs "${fqdir}/*_R2_*.fastq.gz" "${fqdir}/${base}_R2.fastq.gz"

  # remove chunks (comment out if you prefer to keep them)
  find "${fqdir}" -maxdepth 1 -type f -name "*_I1_*.fastq.gz" -delete 2>/dev/null || true
  find "${fqdir}" -maxdepth 1 -type f -name "*_R1_*.fastq.gz" -delete
  find "${fqdir}" -maxdepth 1 -type f -name "*_R2_*.fastq.gz" -delete

  # sanity check
  if [[ -f "${fqdir}/${base}_R1.fastq.gz" && -f "${fqdir}/${base}_R2.fastq.gz" ]]; then
    local r1 r2
    r1=$(zcat "${fqdir}/${base}_R1.fastq.gz" | wc -l)
    r2=$(zcat "${fqdir}/${base}_R2.fastq.gz" | wc -l)
    if [[ $((r1%4)) -ne 0 || $((r2%4)) -ne 0 || "$r1" -ne "$r2" ]]; then
      echo "[WARN] Read count mismatch for ${base}: R1=${r1}, R2=${r2} (lines)"
    else
      echo "[OK] ${base} FASTQs look consistent"
    fi
  fi
}

starsolo_one_sample() {
  local base="$1"
  local fqdir="${FASTQ_DIR}/${base}"
  local outprefix="${STAR_OUT}/${base}_"

  # resume-safe: if STAR already finished, skip
  if [[ -f "${outprefix}Log.final.out" ]]; then
    echo "[SKIP] STARsolo already done for ${base}"
    return 0
  fi

  echo "[STAR] ${base}"
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

run_limited() {  # run_limited <max_jobs> <cmd> [args...]
  local max="$1"; shift
  while (( $(jobs -rp | wc -l) >= max )); do sleep 1; done
  "$@" &
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

  # 1) BAM → FASTQ (parallel, resume-safe)
  for bam in "${bams[@]}"; do
    run_limited "${CONVERT_JOBS}" convert_one_bam "${bam}"
  done
  wait

  # 2) STARsolo (parallel, resume-safe)
  for bam in "${bams[@]}"; do
    base="$(basename "${bam}" .bam)"
    run_limited "${ALIGN_JOBS}" starsolo_one_sample "${base}"
  done
  wait

  echo "[DONE] All samples processed."
}

main "$@"
