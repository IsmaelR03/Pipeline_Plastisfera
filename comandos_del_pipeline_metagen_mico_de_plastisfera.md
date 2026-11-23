# Guía Paso a Paso del Pipeline Bioinformático para Caracterización Metagenómica de la Plastisfera

Este documento describe de forma detallada y reproducible la instalación, verificación y ejecución del flujo de trabajo bioinformático diseñado para analizar lecturas metagenómicas asociadas a la plastisfera (desde el control de calidad hasta la anotación funcional de enzimas biodegradadoras y genes de resistencia antimicrobiana).

## Fase 0: Preparación del Entorno, Instalación y Comprobación de Programas

Para garantizar la reproducibilidad y evitar conflictos con las librerías del sistema operativo (Linux/WSL Ubuntu), todas las herramientas se instalan dentro de un entorno aislado de **Conda / Bioconda**.

### 1. Actualización de certificados SSL del sistema (WSL/Ubuntu)

Evita fallos de verificación de certificados (`X509`) durante las descargas HTTPS desde repositorios públicos como NCBI/ENA:

```
sudo apt-get update && sudo apt-get install -y ca-certificates
sudo update-ca-certificates

```

### 2. Creación y activación del entorno Conda (`ensamblaje`)

Si aún no existe el entorno, se crea e inicia:

```
conda create -n ensamblaje -y
conda activate ensamblaje

```

### 3. Instalación de las herramientas del pipeline

Se instalan en el entorno activo las herramientas de descarga (`sra-tools`, `entrez-direct`), control de calidad (`fastqc`), filtrado (`fastp`), predicción génica (`prodigal`) y alineamiento por homología (`blast`):

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

### 4. Comprobación de que todos los programas están instalados correctamente

Antes de iniciar el procesamiento de datos, se verifica la disponibilidad y versión de cada ejecutable en el entorno:

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

## Fase 1: Adquisición y Extracción de Secuencias Metagenómicas

### 1. Creación de la carpeta de trabajo

```
mkdir -p ~/plastisfera_datoos
cd ~/plastisfera_datoos

```

### 2. Descarga desde NCBI SRA/DDBJ (`prefetch`) y extracción a FASTQ (`fasterq-dump`)

Se descarga el archivo `.sra` de la accesión seleccionada (por ejemplo, el experimento `DRX944086`, que resuelve la corrida `DRR966448`), se extraen las lecturas en formato `.fastq` y se comprimen con `gzip` para optimizar almacenamiento y compatibilidad:

```
prefetch DRX944086
fasterq-dump --split-files --threads 4 --progress DRR966448/DRR966448.sra -O .
gzip -f DRR966448*.fastq
ls -lh DRR966448*.fastq.gz

```

* **Nota:** El parámetro `--split-files` separa automáticamente las lecturas en `_1.fastq` y `_2.fastq` si la librería es *paired-end*, o genera `_1.fastq` si es *single-end*.

## Fase 2: Estructura de Directorios del Análisis

Se organiza el proyecto en subdirectorios numerados para mantener trazabilidad en cada etapa:

```
mkdir -p prueba_pipeline/{01_qc_raw,02_clean,03_qc_clean,04_orfs,05_referencias,06_blastp}

```

## Fase 3: Control de Calidad Inicial (`FastQC`)

Evalúa la calidad por base (puntajes Phred), contenido GC, presencia de adaptadores y secuencias sobrerrepresentadas en las lecturas crudas:

```
fastqc -t 4 --outdir prueba_pipeline/01_qc_raw DRR966448*.fastq.gz

```

* **`-t 4`**: Utiliza 4 hilos de procesamiento.

* **`--outdir`**: Especifica el directorio donde se guardan los reportes `.html` y `.zip`.

## Fase 4: Filtrado y Recorte Automatizado (`fastp`)

Elimina lecturas de baja calidad, adaptadores residuales y fragmentos demasiado cortos (`< 50` pb), generando archivos depurados y reportes cuantitativos en JSON y HTML. El bloque condicional detecta automáticamente si la muestra es *paired-end* o *single-end*:

```
if [ -f "DRR966448_2.fastq.gz" ]; then
  fastp \
    -i DRR966448_1.fastq.gz -I DRR966448_2.fastq.gz \
    -o prueba_pipeline/02_clean/DRR966448_clean_1.fastq.gz \
    -O prueba_pipeline/02_clean/DRR966448_clean_2.fastq.gz \
    -j prueba_pipeline/02_clean/fastp_DRR966448.json \
    -h prueba_pipeline/02_clean/fastp_DRR966448.html \
    --detect_adapter_for_pe --qualified_quality_phred 20 --length_required 50 --thread 4
  LIMPIA="prueba_pipeline/02_clean/DRR966448_clean_1.fastq.gz"
else
  fastp \
    -i DRR966448_1.fastq.gz \
    -o prueba_pipeline/02_clean/DRR966448_clean_1.fastq.gz \
    -j prueba_pipeline/02_clean/fastp_DRR966448.json \
    -h prueba_pipeline/02_clean/fastp_DRR966448.html \
    --qualified_quality_phred 20 --length_required 50 --thread 4
  LIMPIA="prueba_pipeline/02_clean/DRR966448_clean_1.fastq.gz"
fi

```

* **`--qualified_quality_phred 20`**: Exige un puntaje de calidad Phred $\ge 20$.

* **`--length_required 50`**: Descarta lecturas con longitud inferior a 50 pares de bases.

## Fase 5: Control de Calidad Post-Filtrado (`FastQC`)

Verifica que las lecturas depuradas cumplan los criterios de calidad antes de pasar a la predicción de genes:

```
fastqc -t 4 --outdir prueba_pipeline/03_qc_clean prueba_pipeline/02_clean/DRR966448_clean*.fastq.gz

```

## Fase 6: Predicción de Marcos Abiertos de Lectura (`Prodigal`)

Convierte las lecturas limpias de formato FASTQ a FASTA y ejecuta `Prodigal` en modo metagenómico (`-p meta`) para predecir secuencias codificantes (ORFs) y traducirlas a aminoácidos (`.faa`):

```
gunzip -c "${LIMPIA}" | awk 'NR%4==1 {print ">" substr($0,2)} NR%4==2 {print}' > prueba_pipeline/04_orfs/lecturas_limpias.fasta

prodigal \
  -i prueba_pipeline/04_orfs/lecturas_limpias.fasta \
  -a prueba_pipeline/04_orfs/orfs.faa \
  -d prueba_pipeline/04_orfs/orfs.fna \
  -p meta -q

rm -f prueba_pipeline/04_orfs/lecturas_limpias.fasta

```

* **`-a`**: Archivo de salida con las secuencias proteicas predichas (`.faa`).

* **`-d`**: Archivo de salida con las secuencias nucleotídicas predichas (`.fna`).

* **`-p meta`**: Modo optimizado para comunidades metagenómicas mixtas.

## Fase 7: Base de Referencia y Anotación Funcional (`makeblastdb` + `BLASTP`)

### 1. Descarga acotada de proteínas de referencia desde NCBI Protein

Se realiza una búsqueda acotada mediante `esearch` y `efetch` enfocada en enzimas degradadoras de plásticos (`PETase`, `cutinase`, `alkB`) y genes de resistencia (`sul1`) filtrados por RefSeq (`srcdb_refseq[Properties]`) para evitar cortes por descargas masivas:

```
esearch -db protein -query '(PETase[Title] OR cutinase[Title] OR alkB[Gene] OR sul1[Gene]) AND srcdb_refseq[Properties]' | \
efetch -format fasta > prueba_pipeline/05_referencias/ref_prueba.faa

```

### 2. Indexación de la base de datos BLAST y alineamiento con `BLASTP`

Se formatea el archivo FASTA de referencia en una base de datos proteica local y se alinean los ORFs predichos con un umbral de significancia `e-value <= 1e-5`:

```
makeblastdb \
  -in prueba_pipeline/05_referencias/ref_prueba.faa \
  -dbtype prot \
  -out prueba_pipeline/05_referencias/db_ref_prueba

blastp \
  -query prueba_pipeline/04_orfs/orfs.faa \
  -db prueba_pipeline/05_referencias/db_ref_prueba \
  -out prueba_pipeline/06_blastp/DRR966448_hits.tsv \
  -evalue 1e-5 -num_threads 4 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"

```

## Fase 8: Archivos Generados para los Análisis Posteriores

Al finalizar la ejecución completa del pipeline, la estructura de archivos resultante contiene todos los insumos y reportes listos para su interpretación bioinformática y estadístico-descriptiva:

```
~/plastisfera_datoos/
├── DRR966448/
│   └── DRR966448.sra                              # Archivo crudo comprimido descargado desde SRA/DDBJ
├── DRR966448_1.fastq.gz                           # Lecturas crudas extraídas en formato FASTQ comprimido
└── prueba_pipeline/
    ├── 01_qc_raw/
    │   ├── DRR966448_1_fastqc.html                # Reporte visual de calidad inicial (FastQC)
    │   └── DRR966448_1_fastqc.zip                 # Métricas detalladas de calidad inicial
    ├── 02_clean/
    │   ├── DRR966448_clean_1.fastq.gz             # Lecturas depuradas de alta calidad (Q >= 20, >= 50 pb)
    │   ├── fastp_DRR966448.html                   # Reporte visual comparativo antes/después del filtrado
    │   └── fastp_DRR966448.json                   # Métricas cuantitativas estructuradas de retención y calidad
    ├── 03_qc_clean/
    │   ├── DRR966448_clean_1_fastqc.html          # Reporte visual de calidad post-filtrado (FastQC)
    │   └── DRR966448_clean_1_fastqc.zip           # Métricas detalladas post-filtrado
    ├── 04_orfs/
    │   ├── orfs.faa                               # Secuencias de aminoácidos de los ORFs predichos (Prodigal)
    │   └── orfs.fna                               # Secuencias de nucleótidos de los ORFs predichos (Prodigal)
    ├── 05_referencias/
    │   ├── ref_prueba.faa                         # Secuencias proteicas de referencia descargadas (RefSeq)
    │   └── db_ref_prueba.*                        # Archivos binarios indexados de la base BLAST (.phr, .pin, .psq, etc.)
    └── 06_blastp/
        └── DRR966448_hits.tsv                     # Tabla final de anotación funcional (query, sujeto, % identidad, longitud, e-value, bitscore y descripción)

```

### Inspección rápida de los archivos generados

```
grep -c "^>" prueba_pipeline/04_orfs/orfs.faa
grep -c "^>" prueba_pipeline/05_referencias/ref_prueba.faa
wc -l prueba_pipeline/06_blastp/DRR966448_hits.tsv
head -n 10 prueba_pipeline/06_blastp/DRR966448_hits.tsv

```