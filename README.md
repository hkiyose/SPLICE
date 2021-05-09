# SPLICE

## Requirement

* python3
* minimap2 (v2.17 or higher)

## Input file

FASTQ file

## Usage

### Step 1: Preparation of the reference transcriptome file

Download the reference transcriptome file from the UCSC genome browser (https://genome.ucsc.edu/)
1. Select Table Browser from the Tools tab.
2. Select "Genes and Gene Predictions" in the group section.
3. Select the desired database in the Track section and download the file.

Numbering the exons of the reference transcriptome.  
If you use files from two databases (e.g. GENCODE and RefSeq), sort them by gene name and transcript name, and merge identical transcripts.
```
$ cd <path to SPLICE>
$ sh ref_exonnum.sh <path to reference transcriptome> <path to reference genome sequence (FASTA)> <output directory> 
```

Edit the `configure` file  
* Apply the reference genome file (FASTA) to `REF_GENOME_FA`  
* Apply the output `.exonnum` file to `REF_TRANSCRIPT`   
* Apply the output `.exonnum.fa` file to`REF_TRANSCRIPT_FA`

### Step 2: Annotation to reference transcriptome

```
$ sh SPLICE_annot.sh <path to FASTQ> <output directory>
```

### Step 3: Analysis of expression levels

```
$ sh SPLICE_exp.sh <output directory of Step2> <output directory>
```

## Output

`expression.tsv`

| gene | transcript | known/novel | coding/non-coding | transcript length | novel information | sample1 | sample2 | sampleN |
| :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: |
| GAPDH | NM_001256799.2 | known | coding | full | - | 71 | 50 | 31 |
| CSDE1 | NM_001242891.1 | known | coding | partial | - | 346 | 40 | 88 |
| SAP18 | - | novel_exon_length | coding | - | \*,6,8,10,\*/\*,k,l,k,\*/21140681,21147186/\*,0,0,0,-35,0,0,0,0,\* | 8 | 0 | 9 |

Notation of novel transcript
![Notation of novel transcript](https://github.com/hkiyose/SPLICE/images/novelpng)

## Example

Preparation of the reference transcriptome file
```
$ git clone https://github.com/hkiyose/SPLICE
$ cd SPLICE
$ sh ref_exonnum.sh <path to reference transcriptome file> <path to reference genome sequence (FASTA)> ./example/ref
```

Edit the `configure` file as described above.

Annotation to reference transcriptome.
```
$ sh SPLICE_annot.sh ./example/fastq/sample1_test.fastq ./example/annot
$ sh SPLICE_annot.sh ./example/fastq/sample2_test.fastq ./example/annot
```

Analysis of expression levels
```
$ sh SPLICE_exp.sh ./example/annot ./example/exp
```

## Parameter settings in configuration file
If you want to use different parameters, please change the configuration file. 

`BQ_FILT` - Read quality cutoff. Minimum average base quality score (15)  
`MIN_SC_LEN` - Minimum length of the softclip region to be remapped (60)  
`MQ_FILT` - Mapping quality cutoff (0)  
`MIN_FUSION_DIST` - Minimum distance of each transcript in the fusion transcript (200000)  
`MIN_FUSION_READ` - Minimum number of support reads for fusion transcripts (1)  
`MIN_FUSION_FREQ` - Minimum frequency of the fusion transcript in the total amount of the gene (%) (0.1)  
`MQ_FILT_NOVEL_EXON` - Mapping quality cutoff of novel exon (1)  
`MIN_READ_NUM` - Minimum number of support reads (3)  
`MIN_READ_FREQ` - Minimum frequency of the transcript in the total amount of the gene (%) (1)  
`RANGE_SJ_EVA` - Maximum change of novel exon length for evaluate the error rate (20)
`ERR_RATE_FILT` - Mapping error rate cutoff at splicing junctino sites (%) (20)  
`MIN_NOVEL_LEN_GAP` - Minimum change of novel exon length (5)  
`MIN_NOVEL_EXON_LEN` - Minimum length of novel exon (60)  

## License
GPLv3

## Contact
Hiroki Kiyose - hkiyose@m.u-tokyo.ac.jp


