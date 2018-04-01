# CAFE

Composition Anomaly and Feature Enrichment (CAFE), is genomic island prediction tool that utilizes sequence composition and functional information to identify genomic islands.

sudo chmod 775 cafe
sudo chmod 775 cafe.out

## Example
First set the permissions for the file
```
sudo chmod 775 cafe
sudo chmod 775 cafe.out
```
To run program on example files
```
./cafe -phylo -genus bartonella -gbk example.gbk
```

## Input
CAFE requires genome sequence and annotations to predict genomic islands. These can be provided as a single Genbank file.
```
./cafe.pl [options] -gbk  genome.gbk
```
If the genome is not annotated then CAFE can identify marker genes for genomic island. This requires Prodigal and Hmmer be installed and included in path
```
./cafe.pl [options] -annot genome.fna
```
To use the Phylogenetic module of CAFE, first download reference protein sequence files (.faa format) in faa folder. CAFE requires atleast reference protein files for comparison. faa folder should only have reference protein sequence files, remove pre existing files in faa folder if not running on example genome. 
```
./cafe.pl [options] -phylo -genus [genus_name] -gbk genome.gbk
```
Note this requires BLAST version 2.6 or higher. Phylogenetic module also requires that users specify the name 
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





