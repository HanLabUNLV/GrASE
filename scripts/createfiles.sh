!/bin/bash
set -eux

outdir=${HOME}/GrASE_simulation

mkdir -p ${outdir}/gtf
# Split gtf files by gene
awk -v outdir="${outdir}/gtf" '{
    if (match($0, /gene_id "[^"]+"/)) {
        id = substr($0, RSTART+9, RLENGTH-10);
        f = outdir "/" id ".gtf";
        print >> f;
        close(f);
    }
}' ${outdir}/ref/gencode.v28.annotation.gtf

mkdir -p ${outdir}/ref
# Create genelist from the files we just made 
ls ${outdir}/gtf | sed 's/\.gtf$//' > ${outdir}/ref/genelist

mkdir -p ${outdir}/dexseq.gff
# Proceed with dexseq prep for each gene
cat ${outdir}/ref/genelist | parallel -j 30 -I {} "python dexseq_prepare_annotation.py ${outdir}/gtf/{}.gtf ${outdir}/dexseq.gff/{}.dexseq.gff"


# concat dexseq gff into one file (this is to generate true nonaggregated dexseq gff without skipped exons)
find ${outdir}/dexseq.gff/ -name "*.dexseq.gff" -exec cat {} + > ${outdir}/ref/gencode.v28.dexseq.bygene.cat.gff
sort -k1,1 -k4,4n ${outdir}/ref/gencode.v28.dexseq.bygene.cat.gff > ${outdir}/ref/gencode.v28.dexseq.bygene.gff
rm ${outdir}/ref/gencode.v28.dexseq.bygene.cat.gff

# Now detect overlapping genes and determine which to filter out (prioritizing protein-coding and more exonic parts)
# we will still run the analysis on the overlapping genes, but take into account the overlaps in the interpretation
python3 filter_overlapping_genes.py  ~/GrASE_simulation/ref/gencode.v28.dexseq.bygene.gff  ~/GrASE_simulation/gtf  overlapping_genes &> ~/GrASE_simulation/ref/overlap.summary.log

# generate GrASE graphs for each gene 
Rscript ../Rpkg/scripts/generate_graphs.R --indir ~/GrASE_simulation 


