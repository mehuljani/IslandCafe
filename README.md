# CAFE

Composition Anomaly and Feature Enrichment (CAFE), is genomic island prediction tool that utilizes sequence composition and functional information to identify genomic islands.


## Input
CAFE requires genome sequence and annotations to predict genomic islands. These can be provided in one of the following ways:
*Fasta file containing the genome sequence and another file containing annotations in ptt format.
```
./cafe.pl [options] genome.fna annotation.ptt
```
*A single file containing both the sequence and annotations in gbk format.
```
./cafe.pl [options] -gbk  genome.gbk
```
*If the genome is not annotated then CAFE can identify marker genes for genomic island. This requires Prodigal and Hmmer be installed and included in path
```
./cafe.pl [options] -annot genome.fna
```

### Optional Input
-annot    Annotate marker genes (This requires Prodigal and Hmmer are installed and in path)
-Thres    Provide segmentation, contiguous clustering and non-contiguous clustering thresholds (range: 0-1)
-gbk      Use genbank as input file
-out      Output file name
-verbose  print on screen
-expert   keep temporary files for user analyses
-visual   Make a map of genomic islands (Requires CGView be installed and in path)


## Output
CAFE outputs two tab separated text files. File with suffix CAFE_full_version.txt shows genomic island predictions using CAFE full version. Likewise, file with suffix CAFE_marker_version.txt. shows genomic island predictions using CAFE marker version.
Both files has four columns showing genomic island id, start and end co-ordinates, and length of the genomic island.

## Example
./cafe.pl -out outfile -verbose -gbk example.gbk




