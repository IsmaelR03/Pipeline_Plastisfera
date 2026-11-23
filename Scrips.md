```
sudo apt-get update && sudo apt-get install -y ca-certificates
sudo update-ca-certificates
```
```
conda create -n ensamblaje -y
conda activate ensamblaje
```
```
conda install -c bioconda -c conda-forge \
  sra-tools \
  fastqc \
  fastp \
  prodigal \
  blast \
  entrez-direct \
  ca-certificates \
  openssl -y
```
```
prefetch --version
fasterq-dump --version
fastqc --version
fastp --version
prodigal -v
makeblastdb -version
blastp -version
esearch -version
```
```
mkdir -p ~/plastisfera_datos
cd ~/plastisfera_datos
prefetch DRX944086
fasterq-dump --split-files --threads 4 --progress DRR966448/DRR966448.sra -O .
gzip -f DRR966448*.fastq
ls -lh DRR966448*.fastq.gz
```
```
mkdir -p Propuesta_pipeline/{01_qc_raw,02_clean,03_qc_clean,04_orfs,05_referencias,06_blastp}
fastqc -t 4 --outdir Propuesta_pipeline/01_qc_raw DRR966448*.fastq.gz
if [ -f "DRR966448_2.fastq.gz" ]; then
  fastp \
    -i DRR966448_1.fastq.gz -I DRR966448_2.fastq.gz \
    -o Propuesta_pipeline/02_clean/DRR966448_clean_1.fastq.gz \
    -O Propuesta_pipeline/02_clean/DRR966448_clean_2.fastq.gz \
    -j Propuesta_pipeline/02_clean/fastp_DRR966448.json \
    -h Propuesta_pipeline/02_clean/fastp_DRR966448.html \
    --detect_adapter_for_pe --qualified_quality_phred 20 --length_required 50 --thread 4
  LIMPIA="Propuesta_pipeline/02_clean/DRR966448_clean_1.fastq.gz"
else
  fastp \
    -i DRR966448_1.fastq.gz \
    -o Propuesta_pipeline/02_clean/DRR966448_clean_1.fastq.gz \
    -j Propuesta_pipeline/02_clean/fastp_DRR966448.json \
    -h Propuesta_pipeline/02_clean/fastp_DRR966448.html \
    --qualified_quality_phred 20 --length_required 50 --thread 4
  LIMPIA="Propuesta_pipeline/02_clean/DRR966448_clean_1.fastq.gz"
fi
```
```
fastqc -t 4 --outdir Propuesta_pipeline/03_qc_clean Propuesta_pipeline/02_clean/DRR966448_clean*.fastq.gz
gunzip -c "${LIMPIA}" | awk 'NR%4==1 {print ">" substr($0,2)} NR%4==2 {print}' > Propuesta_pipeline/04_orfs/lecturas_limpias.fasta
prodigal \
  -i Propuesta_pipeline/04_orfs/lecturas_limpias.fasta \
  -a Propuesta_pipeline/04_orfs/orfs.faa \
  -d Propuesta_pipeline/04_orfs/orfs.fna \
  -p meta -q
```
```
esearch -db protein -query '(PETase[Title] OR cutinase[Title] OR alkB[Gene] OR sul1[Gene]) AND srcdb_refseq[Properties]' | \
efetch -format fasta > Propuesta_pipeline/05_referencias/ref_prueba.faa
```
```
makeblastdb \
  -in Propuesta_pipeline/05_referencias/ref_prueba.faa \
  -dbtype prot \
  -out Propuesta_pipeline/05_referencias/db_ref_prueba
```
```
blastp \
  -query Propuesta_pipeline/04_orfs/orfs.faa \
  -db Propuesta_pipeline/05_referencias/db_ref_prueba \
  -out Propuesta_pipeline/06_blastp/DRR966448_hits.tsv \
  -evalue 1e-5 -num_threads 4 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
```
```
grep -c "^>" Propuesta_pipeline/04_orfs/orfs.faa
grep -c "^>" Propuesta_pipeline/05_referencias/ref_prueba.faa
wc -l Propuesta_pipeline/06_blastp/DRR966448_hits.tsv
head -n 10 Propuesta_pipeline/06_blastp/DRR966448_hits.tsv
wc -l Propuesta_pipeline/06_blastp/DRR966448_hits.tsv
head -n 10 Propuesta_pipeline/06_blastp/DRR966448_hits.tsv
```
