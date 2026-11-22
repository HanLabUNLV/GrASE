#!/bin/bash
# Two-pass STAR alignment: index, pass 1 (junction discovery), junction filter,
# index with the filtered junctions, pass 2 (BAMs). Each step skips work already
# done, so a stopped run resumes where it left off.
#
# Usage:
#   bash star_two_pass.sh --samples FILE --fasta FILE --gtf FILE --outdir DIR \
#     --sjdb_overhang N --overhang_min N [options]
#
# Required:
#   --samples        TSV: sample, group, fastq1[, fastq2]; one line per sample.
#                    Junctions are filtered per group, then merged across groups.
#   --fasta          genome FASTA (primary assembly: no alternate haplotypes or
#                    patches, which make reads from those loci multimap)
#   --gtf            annotation for both indexes
#   --outdir         output root
#   --sjdb_overhang  read length - 1
#   --overhang_min   --alignSJDBoverhangMin of both passes
#
# Options:
#   --star           STAR binary [STAR]
#   --samtools       samtools, to index the pass 2 BAMs (MAJIQ needs the index) [samtools]
#   --threads        threads per alignment [4]
#   --tmpdir         directory for STAR's per-sample temporary files (--outTmpDir), e.g.
#                    on local NVMe; removed after each sample [<sample's output dir>/tmp_<sample>]
#   --jobs           alignments at once [12]; the genome is loaded once into
#                    shared memory and used by all of them
#   --index_threads  threads for genomeGenerate [24]
#   --steps          comma list of index1,pass1,sj,index2,pass2 [all]
#   --min_samples, --min_unique, --min_maxoverhang, --final_unique
#                    junction filter per group [5, 10, 10, 20]: canonical motif,
#                    found in >= min_samples samples with >= min_unique unique
#                    reads in total and a maximum overhang >= min_maxoverhang;
#                    then >= final_unique unique reads
#
# Output, under --outdir:
#   index_pass1/, pass1/<group>/<sample>.SJ.out.tab, sj/<group>.{raw,highconf,final}.SJ.tab,
#   sj/all.final.SJ.tab, index_pass2/,
#   bam/<group>/pass2_<sample>.{Aligned.sortedByCoord.out.bam(.bai),SJ.out.tab,Log.final.out}
#   run_info.txt: STAR version, inputs and their md5sums
set -euo pipefail

STAR=STAR; SAMTOOLS=samtools; TMPDIR_=; THREADS=4; JOBS=12; ITHREADS=24; STEPS=index1,pass1,sj,index2,pass2
MIN_SAMPLES=5; MIN_UNIQUE=10; MIN_MAXOH=10; FINAL_UNIQUE=20
while [ $# -gt 0 ]; do
  case "$1" in
    --samples) SAMPLES=$2 ;; --fasta) FASTA=$2 ;; --gtf) GTF=$2 ;; --outdir) OUT=$2 ;;
    --sjdb_overhang) SJDB_OH=$2 ;; --overhang_min) OH_MIN=$2 ;; --star) STAR=$2 ;; --samtools) SAMTOOLS=$2 ;;
    --threads) THREADS=$2 ;; --tmpdir) TMPDIR_=$2 ;; --jobs) JOBS=$2 ;; --index_threads) ITHREADS=$2 ;; --steps) STEPS=$2 ;;
    --min_samples) MIN_SAMPLES=$2 ;; --min_unique) MIN_UNIQUE=$2 ;;
    --min_maxoverhang) MIN_MAXOH=$2 ;; --final_unique) FINAL_UNIQUE=$2 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
  shift 2
done
for v in SAMPLES FASTA GTF OUT SJDB_OH OH_MIN; do
  [ -n "${!v:-}" ] || { echo "missing a required option (see the header of $0)" >&2; exit 1; }
done
for f in "$SAMPLES" "$FASTA" "$GTF"; do [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }; done
"$STAR" --version > /dev/null || { echo "$STAR does not run" >&2; exit 1; }
"$SAMTOOLS" --version > /dev/null || { echo "$SAMTOOLS does not run" >&2; exit 1; }
mkdir -p "$OUT"; OUT=$(realpath "$OUT")
[ -z "$TMPDIR_" ] || { mkdir -p "$TMPDIR_" && TMPDIR_=$(realpath "$TMPDIR_"); } || { echo "cannot create --tmpdir $TMPDIR_" >&2; exit 1; }
has () { [[ ",$STEPS," == *",$1,"* ]]; }
I1=$OUT/index_pass1; I2=$OUT/index_pass2
SJALL=$OUT/sj/all.final.SJ.tab

{ echo "STAR $("$STAR" --version)"; echo "fasta $FASTA $(md5sum < "$FASTA" | cut -d' ' -f1)"
  echo "gtf $GTF $(md5sum < "$GTF" | cut -d' ' -f1)"; echo "samples $SAMPLES $(wc -l < "$SAMPLES")"
  echo "sjdbOverhang $SJDB_OH alignSJDBoverhangMin $OH_MIN"
  echo "junction filter: samples >= $MIN_SAMPLES, unique >= $MIN_UNIQUE, max overhang >= $MIN_MAXOH, final unique >= $FINAL_UNIQUE"
} > "$OUT/run_info.txt"

index () {  # index <dir> [junction file]
  local d=$1 sj=${2:-}
  [ -s "$d/SAindex" ] && { echo "  $d exists"; return 0; }
  mkdir -p "$d"
  "$STAR" --runThreadN "$ITHREADS" --runMode genomeGenerate --genomeDir "$d" \
    --genomeFastaFiles "$FASTA" --sjdbGTFfile "$GTF" --sjdbOverhang "$SJDB_OH" \
    ${sj:+--sjdbFileChrStartEnd "$sj"} --outFileNamePrefix "$d/"
}

# Alignment settings of both passes, as in the DICE and simulation runs
ALIGN=(--alignIntronMax 299999 --alignSJDBoverhangMin "$OH_MIN" --alignEndsType EndToEnd
       --outFilterMismatchNmax 3 --chimSegmentMin 2)

align_one () {  # align_one <pass> <index> <sample> <group> <fastq1> [fastq2]
  local pass=$1 idx=$2 s=$3 g=$4 f1=$5 f2=${6:-} d pre tmp rf=()
  if [ "$pass" = 1 ]; then d=$OUT/pass1/$g; pre=$d/$s.; else d=$OUT/bam/$g; pre=$d/pass2_$s.; fi
  [ -s "${pre}Log.final.out" ] && return 0
  mkdir -p "$d"
  # STAR wants a tmp dir that does not exist yet; one left by an interrupted run is removed
  if [ -n "$TMPDIR_" ]; then tmp=$TMPDIR_/pass${pass}_${g}_$s; else tmp=$d/tmp_$s; fi
  rm -rf "$tmp"
  [[ "$f1" == *.gz ]] && rf=(--readFilesCommand zcat)
  local out=(--outSAMtype None)
  [ "$pass" = 2 ] && out=(--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 8000000000 --outSAMstrandField intronMotif)
  "$STAR" --runThreadN "$THREADS" --genomeDir "$idx" --genomeLoad LoadAndKeep "${ALIGN[@]}" \
    --readFilesIn "$f1" ${f2:+"$f2"} "${rf[@]}" "${out[@]}" \
    --outFileNamePrefix "$pre" --outTmpDir "$tmp" > "${pre}stdout" 2>&1 \
    || { rm -f "${pre}Log.final.out"; rm -rf "$tmp"; echo "FAILED $s (pass $pass): see ${pre}Log.out" >&2; return 1; }
  rm -rf "$tmp"
  [ "$pass" = 1 ] && mv "${pre}SJ.out.tab" "$d/$s.SJ.out.tab"
  if [ "$pass" = 2 ]; then
    "$SAMTOOLS" index "${pre}Aligned.sortedByCoord.out.bam" \
      || { rm -f "${pre}Log.final.out"; echo "FAILED $s: samtools index" >&2; return 1; }
  fi
  return 0
}
export -f align_one; export OUT STAR SAMTOOLS THREADS TMPDIR_; export ALIGN_STR="${ALIGN[*]}"

run_pass () {  # run_pass <pass> <index>
  local idx=$2
  "$STAR" --genomeDir "$idx" --genomeLoad LoadAndExit --outFileNamePrefix "$OUT/load_$1." > /dev/null
  trap '"$STAR" --genomeDir "'"$idx"'" --genomeLoad Remove --outFileNamePrefix "'"$OUT"'/unload." > /dev/null 2>&1 || true' EXIT
  # ALIGN is an array, which export cannot carry; rebuild it in each job
  awk -F'\t' 'NF >= 3' "$SAMPLES" | xargs -P "$JOBS" -L 1 bash -c \
    'ALIGN=($ALIGN_STR); align_one '"$1"' '"$idx"' "$@"' _
  "$STAR" --genomeDir "$idx" --genomeLoad Remove --outFileNamePrefix "$OUT/unload_$1." > /dev/null
  trap - EXIT
  local n=$(awk -F'\t' 'NF >= 3' "$SAMPLES" | wc -l) done_=0
  if [ "$1" = 1 ]; then done_=$(find "$OUT/pass1" -name '*.Log.final.out' | wc -l)
  else done_=$(find "$OUT/bam" -name 'pass2_*.Log.final.out' | wc -l); fi
  echo "  pass $1: $done_ of $n samples aligned"
  [ "$done_" -eq "$n" ] || { echo "pass $1 incomplete" >&2; exit 1; }
}

if has index1; then echo "=== index, pass 1  $(date '+%F %T')"; index "$I1"; fi
if has pass1; then echo "=== pass 1  $(date '+%F %T')"; run_pass 1 "$I1"; fi

if has sj; then
  echo "=== junction filter  $(date '+%F %T')"
  mkdir -p "$OUT/sj"
  for g in $(cut -f2 "$SAMPLES" | sort -u); do
    cat "$OUT/pass1/$g"/*.SJ.out.tab | awk '$5 > 0' | sort -k1,1 -k2,2n -k3,3n | uniq > "$OUT/sj/$g.raw.SJ.tab"
    # per junction: samples it occurs in, unique reads summed, largest overhang
    awk -F'\t' -v OFS='\t' -v ms="$MIN_SAMPLES" -v mu="$MIN_UNIQUE" -v mo="$MIN_MAXOH" '
      $5 == 1 || $5 == 2 { k = $1 FS $2 FS $3 FS $4; n[k]++; u[k] += $7; m[k] = $5; a[k] = $6
                           if ($9 > o[k]) o[k] = $9 }
      END { for (k in n) if (n[k] >= ms && u[k] >= mu && o[k] >= mo) { split(k, x, FS)
              print x[1], x[2], x[3], x[4], m[k], a[k], u[k], o[k] } }' "$OUT/sj/$g.raw.SJ.tab" \
      | sort -k1,1 -k2,2n -k3,3n > "$OUT/sj/$g.highconf.SJ.tab"
    awk -F'\t' -v fu="$FINAL_UNIQUE" -v mo="$MIN_MAXOH" '$7 >= fu && $8 >= mo' "$OUT/sj/$g.highconf.SJ.tab" > "$OUT/sj/$g.final.SJ.tab"
    echo "  $g: $(wc -l < "$OUT/sj/$g.final.SJ.tab") junctions"
  done
  cat "$OUT/sj"/*.final.SJ.tab | cut -f1-6 | sort -k1,1 -k2,2n -k3,3n | uniq > "$SJALL"
  echo "  merged: $(wc -l < "$SJALL") junctions"
fi

if has index2; then echo "=== index, pass 2  $(date '+%F %T')"; [ -s "$SJALL" ] || { echo "missing $SJALL" >&2; exit 1; }; index "$I2" "$SJALL"; fi
if has pass2; then echo "=== pass 2  $(date '+%F %T')"; run_pass 2 "$I2"; fi
echo "=== done  $(date '+%F %T')"
