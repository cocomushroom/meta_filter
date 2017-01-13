# Metatranscriptome filtering for nutrient-genes of fungal community
Metatranscriptomic data of symbiotic system are complicated to trace certain group of organisms or genes of interst.
Files in this folder are created to deal with complex metatranscriptomic data. I'm interested in nutrient-related genes in fungal communities associated with the plant host. 

## Majors steps

### Data assembly
- Quality control (trimmomatic)
- Filter and remove rRNA (sortMeRNA)
- Assembly (Trinity) (kept gene level)
- Assembly evaluation ()
- Clustering (CD-HIT)

### 1st sorting organisms by UniProt
- BLAST (BLASTX) against UniProt
- Keep genes with top hit as fungi (BashScript)
- Converting transcripts to proteins (Transdecoder)

### 2nd sorting organims by NCBI NR database
-BLAST (BLASTP) against NCBI NR
-Sorting organisms (MEGAN)

*To load multiple BLAST results into MEGAN, some modification must be made beside just concatenation.
See

### Annotation
-Annotation based on Pfam, UniProt and NCBI (Trinotate)

### Manually exam for obvious erros, filter or extract gene lists
- Filter out unwant genes (e.g. housekeeping genes)
- Extract wanted genes (e.g. transmembrane transporter)


### Expression examination
- Mapping (BOWTIE2, samtool)
- Concatenating matrix (BashScript)
- Differential Gene Expression test (DEseq2)
- Rank, filter and extract annotation info (Rscript)
