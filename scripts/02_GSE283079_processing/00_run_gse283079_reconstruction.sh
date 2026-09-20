#!/usr/bin/env bash
set -euo pipefail

project_dir="$(pwd)"
manifest_csv="$project_dir/结果/阶段0.8/GSE283079_sample_key_audit.csv"
manifest_tsv="$project_dir/数据/中间/Stage2A_trial/GSE283079_run_manifest.tsv"
raw_dir="$project_dir/数据/原始/Stage2A_SRA"
quant_dir="$project_dir/数据/结果/Stage2A_quant"
log_dir="$project_dir/05_日志/Stage2A独立验证"
tmp_fastq="/private/tmp/stage2a_batch_fastq"
tmp_work="/private/tmp/stage2a_batch_temp"
sra_tool="$project_dir/数据/软件/Stage2A/mamba_root/envs/stage2a/bin"
index_dir="$project_dir/数据/参考/Stage2A_GENCODE47/salmon_index_v47"

mkdir -p "$quant_dir" "$log_dir" "$tmp_fastq" "$tmp_work"
awk -F, 'NR>1 {print $1 "\t" $3 "\t" $7}' "$manifest_csv" > "$manifest_tsv"

while IFS=$'\t' read -r gsm group srr; do
  [ -z "$srr" ] && continue
  quant_sample="$quant_dir/$srr"
  if [ -s "$quant_sample/quant.sf" ]; then
    echo "$(date -Iseconds) SKIP $gsm $group $srr quant_exists" | tee -a "$log_dir/batch_pipeline.log"
    continue
  fi
  echo "$(date -Iseconds) START $gsm $group $srr" | tee -a "$log_dir/batch_pipeline.log"
  sra_file="$raw_dir/$srr/$srr.sra"
  if [ ! -s "$sra_file" ]; then
    "$sra_tool/prefetch" "$srr" --max-size 6G -O "$raw_dir" --progress >> "$log_dir/${srr}_prefetch.log" 2>&1
  fi
  rm -rf "$tmp_fastq/$srr" "$tmp_work/$srr"
  mkdir -p "$tmp_fastq/$srr" "$tmp_work/$srr"
  "$sra_tool/fasterq-dump" --split-files --threads 8 --temp "$tmp_work/$srr" --progress -O "$tmp_fastq/$srr" "$sra_file" > "$log_dir/${srr}_fasterq.log" 2>&1
  "$sra_tool/fastp" -i "$tmp_fastq/$srr/${srr}_1.fastq" -I "$tmp_fastq/$srr/${srr}_2.fastq" -o "$tmp_work/$srr/${srr}_R1.clean.fastq.gz" -O "$tmp_work/$srr/${srr}_R2.clean.fastq.gz" --detect_adapter_for_pe --thread 8 --json "$quant_dir/${srr}_fastp.json" --html "$quant_dir/${srr}_fastp.html" > "$log_dir/${srr}_fastp.log" 2>&1
  "$sra_tool/salmon" quant -i "$index_dir" -l A -1 "$tmp_work/$srr/${srr}_R1.clean.fastq.gz" -2 "$tmp_work/$srr/${srr}_R2.clean.fastq.gz" --seqBias --gcBias --numBootstraps 0 -p 8 -o "$quant_sample" > "$log_dir/${srr}_salmon.log" 2>&1
  rm -rf "$tmp_fastq/$srr" "$tmp_work/$srr"
  echo "$(date -Iseconds) DONE $gsm $group $srr" | tee -a "$log_dir/batch_pipeline.log"
done < "$manifest_tsv"
