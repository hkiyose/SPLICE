# SPLICE

Analysis pipeline for long-read RNA-seq data using Nanopore technology

## Requirement

* python3 (v3.5 or higher)
* minimap2 (v2.17 or higher)
* luigi

## Input file

FASTQ file

## Usage

### Step 1: Preparation of the reference transcriptome file

Download the reference transcriptome file from the UCSC genome browser (https://genome.ucsc.edu/)
1. Select Table Browser from the Tools tab.
2. Select "Genes and Gene Predictions" in the group section.
3. Select the desired database in the Track section and download the file.

Numbering the exons of the reference transcriptome.  
If you use files from two databases (e.g. GENCODE and RefSeq), sort them by gene name and transcript name, and remove redundant transcripts.
```
$ cd <path to SPLICE>
$ python ref_exonnum.py -h
usage: ref_exonnum.py [-h] -i INPUT -g GENOME [-o OUTPUT] [-p PREFIX] [-w WORKERS] [--tmp TMP]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        reference transcriptome (UCSC table browser output file)
  -g GENOME, --genome GENOME
                        reference genome sequence (FASTA)
  -o OUTPUT, --output OUTPUT
                        output directory
  -p PREFIX, --prefix PREFIX
                        prefix for the output files
  -w WORKERS, --workers WORKERS
                        number of threads
  --tmp TMP             temporary directory
```

### Step 2: Annotation to reference transcriptome

```
$ python annot.py -h
usage: annot.py [-h] -i INPUT -g GENOME -r REF -p PREFIX [-o OUTPUT] [--tmp TMP] [-w WORKERS] [--bq_filt BQ_FILT]
                [--min_sc_len MIN_SC_LEN] [--mq_filt MQ_FILT] [--min_fusion_dist MIN_FUSION_DIST]
                [--max_fusion_bp_merge MAX_FUSION_BP_MERGE] [--min_fusion_read MIN_FUSION_READ]
                [--min_fusion_freq MIN_FUSION_FREQ] [--mq_filt_novel_exon MQ_FILT_NOVEL_EXON]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        query sequence data (FASTQ)
  -g GENOME, --genome GENOME
                        reference genome sequence (FASTA)
  -r REF, --ref REF     reference transcriptome file (output file of ref_exonnum.py)
  -p PREFIX, --prefix PREFIX
                        prefix for the output files
  -o OUTPUT, --output OUTPUT
                        output directory
  --tmp TMP             temporary directory
  -w WORKERS, --workers WORKERS
                        number of threads
  --bq_filt BQ_FILT     read quality cutoff. Minimum average base quality score (15)
  --min_sc_len MIN_SC_LEN
                        minimum length of the softclip region to be remapped (60)
  --mq_filt MQ_FILT     mapping quality cutoff (0)
  --min_fusion_dist MIN_FUSION_DIST
                        minimum distance of each transcript in the fusion transcript (200000)
  --max_fusion_bp_merge MAX_FUSION_BP_MERGE
                        maximam distance to merge fusion gene breakpoints (5)
  --min_fusion_read MIN_FUSION_READ
                        minimum number of support reads for fusion transcripts (1)
  --min_fusion_freq MIN_FUSION_FREQ
                        minimum frequency of the fusion transcript in the total amount of the gene (percentage) (0.1)
  --mq_filt_novel_exon MQ_FILT_NOVEL_EXON
                        mapping quality cutoff of novel exon (1)
```

### Step 3: Analysis of expression levels (Option for multiple analysis)

```
$ python exp.py -h
usage: exp.py [-h] -i INPUT -r REF [-o OUTPUT] [--tmp TMP] [-w WORKERS] [--min_read_num MIN_READ_NUM]
              [--min_read_freq MIN_READ_FREQ] [--range_sj_eva RANGE_SJ_EVA] [--err_rate_filt ERR_RATE_FILT]
              [--min_novel_len_gap MIN_NOVEL_LEN_GAP] [--min_novel_exon_len MIN_NOVEL_EXON_LEN]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input directory (output directory of annot.py)
  -r REF, --ref REF     reference transcriptome file (output file of ref_exonnum.py)
  -o OUTPUT, --output OUTPUT
                        output directory ('./out/exp')
  --tmp TMP             temporary directory ('./tmp/exp_tmp')
  -w WORKERS, --workers WORKERS
                        number of threads
  --min_read_num MIN_READ_NUM
                        minimum number of support reads (3)
  --min_read_freq MIN_READ_FREQ
                        minimum frequency of the transcript in the total amount of the gene (percentage) (1)
  --range_sj_eva RANGE_SJ_EVA
                        maximum change of novel exon length for evaluate the error rate (20)
  --err_rate_filt ERR_RATE_FILT
                        mapping error rate cutoff at splicing junctino sites (percentage) (20)
  --min_novel_len_gap MIN_NOVEL_LEN_GAP
                        minimum change of novel exon length (5)
  --min_novel_exon_len MIN_NOVEL_EXON_LEN
                        minimum length of novel exon (60)
```

## Output

`expression.tsv` - Output file of transcript expression levels (number of supporting reads)

| gene | transcript | known/novel | coding/non-coding | transcript length | novel information | sample1 | sample2 | sampleN |
| :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: |
| GAPDH | NM_001256799.2 | known | coding | full | - | 71 | 50 | 31 |
| CSDE1 | NM_001242891.1 | known | coding | partial | - | 346 | 40 | 88 |
| SAP18 | - | novel_exon_length | coding | - | \*,6,8,10,\*/\*,k,l,k,\*/21140681,21147186/\*,0,0,0,-35,0,0,0,0,\* | 8 | 0 | 9 |

Notation of novel transcript  
![Notation of novel transcript](https://github.com/hkiyose/SPLICE/blob/master/images/novel.png)

`.fusion` - Output file of fusion transcript expression levels (number of supporting reads) for each sample
| Number of reads | Read frequency(%) | GeneA/B | ChrA/B | BreakpointA/B | Read IDs |
| :----: | :----: | :----: | :----: | :----: | :----: |
| 150 | 28.571 | UQCRFS1/YWHAE | chr19/chr17 | 29207585-29207585/1364879-1364879 | 86ada1... |
| 29 | 7.143 | RPS6KA5/TMSB4X | chr14/chrX | 91060446-91060446/12977041-12977041 | 38a936... |
| 12 | 5.128 | TMED10/VPS4B | chr14/chr18 | 75132140-75132140/63390367-63390367 | 1e0cbc... |

breakpoint indicates the range after merging the neighboring breakpoints

## Example

Preparation of the reference transcriptome file
```
$ git clone https://github.com/hkiyose/SPLICE
$ cd SPLICE
$ python ref_exonnum.py -i ./example/gencode_v29_chr1.tsv -g <path to reference genome sequence (FASTA)>
```

Annotation to reference transcriptome.
```
$ python annot.py -i ./example/fastq/sample1_test.fastq -g <path to reference genome sequence (FASTA)> -r ./out/ref_exonnum -p sample1 --tmp ./tmp/annot_tmp/sample1
$ python annot.py -i ./example/fastq/sample2_test.fastq -g <path to reference genome sequence (FASTA)> -r ./out/ref_exonnum -p sample2 --tmp ./tmp/annot_tmp/sample2
```

Analysis of expression levels
```
$ python exp.py -i ./out/annot -r ./out/ref_exonnum
```

## Installation and usage via Docker
Install Docker in your computer, and build a Docker image with the following commands.
```
$ git clone https://github.com/hkiyose/SPLICE.git
$ cd <path to SPLICE>
$ docker build -t splice .
```

The following command mounts the host directory containing the data to the container.
Refer to Step 1 of Usage to download the reference data.
```
$ docker run --rm -it \
  -v <path to directory of reference transcriptome and referece genome sequence (FASTA)>:/ref \
  -v <path to directory of sample data(FASTQ)>:/sample_input \
  -v <path to output directory>:/sample_output \
  splice
```
Then Run according to Usage.

## License
GPLv3

## Contact
Hiroki Kiyose - hkiyose@m.u-tokyo.ac.jp


