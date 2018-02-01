# CAFE

## Required Input
CAFE requires a file containing the genomic sequence in fasta format and annotation file in ptt format
```
./cafe.pl [options] genome.fna annotation.ptt
```
OR

The users can provide a file conatining both, sequence and annotation in genbank (.gbk) format
```
./cafe.pl [options] -annot genome.fna
```

OR

If the annotation file is not available, CAFE will identify marker genes for genomic islands
```
./cafe.pl [options] -annot genome.fna
```
## Optional Input
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




