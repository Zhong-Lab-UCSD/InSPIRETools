# SMOPINseqTools

## Overview
Small molecule-protein interaction sequencing (SMOPIN-seq) is a high-throughput sequencing technology that efficiently detects small molecule-protein associations in vitro. Here, we distribute SMOPINseqTools, a standardized data processing pipeline to identify small molecule-protein associations from fastq files of SMOPIN-seq experiments.<br />

## Workflow
1. Raw read pairs from the PRIM-seq experiment are present in `.fastq` files.
2. The reads are first searched for small molecule library barcode with 0 error tolerance.
3. The reads with library barcode are then searched for 3 cycles of moledule barcodes with 1 base editing tolerance.
4. The read-ends with library and molecule barcodes are assigned as the small-molecule-end reads. The corresponding other ends eassigned as the protein-end reads.
5. Cutadpt is applied to remove 3' linker sequences and 5' adapter sequences from the protein-end reads. 
6. Fastp is then applied to the adapter-trimmed protein-end reads to remove low-quality reads whose mean quality is lower than Q20 and too short reads whose length is shorter than 20 bp.
7. The remaining protein-end reads are mapped to transcriptome with BWA with default parameters. 
8. The mapped protein-end reads are output in `.bed` file with aligned genes and transcriptome alignment information.
6. The protein-end reads are paired with mapped small-molecule reads by read ids. Deduplications are then performed based on UMIs.
12. The kept small-molecule-protein pairs are output as small molecule-protein associations in `SmoProteinAssociations.csv`.


## Software Requirements
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) (2.5 or later)
- [fastp](https://github.com/OpenGene/fastp) (0.22.0 or later)
- [bwa](https://github.com/lh3/bwa) (0.7.17-r1188 or later)
- [samtools](http://www.htslib.org/) (1.6 or later)
- [bedtools](https://bedtools.readthedocs.io/en/latest/) (2.30.0 or later)
- Python 3.4 or later, the following python libraries are required:<br />
    - sys
    - collections
    - cigar
    - glob
    - scipy
    - datetime
## Additional files required
**BWA Index of the transcriptome to be aligned**<br />
You will need to download or build the bwa index of the target trancriptome for PROPERseqTools to use. Here we provide the compressed bwa index built from [RefSeq GRCh38 transcriptome](https://drive.google.com/file/d/1lAV-dVVwVaPi-qVLXibaAvgjtd1e-QMT/view?usp=sharing)

**Transcript, gene and gene type dictionary file**<br />
You will also need a dictionary file that contains the information of transcript ids to their corresponding gene names/gene ids and corresponding gene types in a csv format with the first column being transcript ids, the second column being gene names and the third column being gene types. Here we provide an example dictionary file for [RefSeq GRCh38 genome](https://github.com/Zhong-Lab-UCSD/PRIMseqTools/blob/main/refSeq_tx_gene_type.csv)



## Usage
**Installation**
1. Clone the current github repository to your local machine. For example<br />
`git clone https://github.com/Zhong-Lab-UCSD/PRIMseqTools`
2. Add the following path of the cloned directory to your `.bashrc` file<br />
`export PATH=$PATH:/home/path/to/PRIMseqTools/bin`

**To excute PROPERseqTools, run**
<pre><code>
PRIMseqTools -a /path/to/read1.fastq
             -b /path/to/read2.fastq
             -i /path/to/bwaIndex/transcriptome.fa
             -o /path/to/outputDir
             -g /path/to/refSeq_tx_gene_type.csv
           
</code></pre>



**Required parameters**
<pre><code>
-a     |String, Path to read1 fastq file, fastq.gz also supported
-b     |String, Path to read2 fastq file, fastq.gz also supported
-o     |String, Path to output directory
-i     |String, Path to bwa index of the target transcriptome
-g     |String, Path to transcirpt, gene and gene type dictionary file
</code></pre>
**Other parameters**
<pre><code>
-d     |Float, odds ratio cutoff used to identify RNA-protein associations, default=1
-p     |Float, false discovery rate cutoff used to identify RNA-protein associations, default=0.05
-c     |Float, read count cutoff coefficient used to identify RNA-protein associations, default=4
-j     |String, Job ID to be prepended to the output files and directories, optional, default=PROPERseq"
-t     |Int, Number of working threads, default=2
-r     |Char, (T or F), removal of intermediate files or not, default=T
-h     |Print usage message" 
           
