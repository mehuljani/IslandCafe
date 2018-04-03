# CAFE

Composition Anomaly and Feature Enrichment (CAFE), is genomic island prediction tool that utilizes sequence composition and functional information to identify genomic islands.

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
./cafe [options] -gbk  genome.gbk
```
If the genome is not annotated then CAFE can identify marker genes for genomic island. This requires Prodigal and Hmmer be installed and included in path
```
./cafe [options] -annot genome.fna
```
To use the Phylogenetic module of CAFE, first download reference protein sequence files (.faa format) in faa folder. CAFE requires atleast reference protein files for comparison. faa folder should only have reference protein sequence files, remove pre existing files in faa folder if not running on example genome. Phylogenetic module also requires that genus name be speciied using -genus option
```
./cafe [options] -phylo -genus [genus_name] -gbk genome.gbk
```
Note this requires BLAST version 2.6 or higher. Phylogenetic module also requires that users specify the name 
### Options
--help    Print help and exit

--info    Print program information and exit

--annot   Annotate marker genes. This option is only required if the input file is in fasta format. This option requires prodigal and Hmmer be installed and in path 

--Thres   Provide segmentation, contiguous clustering and non contiguous clustering thresholds (range: 0-1. eg 0.8 0.99999 0.999)

--gbk     Use genbank file as input

--phylo   Use phylogenetic module. Genus must be specified for using phylogenetic module. (eg ./cafe -phylo -genus escherichia -gbk ecoli.gbk)

--genus   Specify genus name of the input genome

--out     Output file name 

--verbose print on screen

--expert  keep temporary files for user analyses

--visual  Make a map of genomic islands (Requires CGView be installed and in path)


## Output
CAFE outputs a tab separated text file. File with suffix CAFE.txt shows genomic island predictions. It has four columns showing genomic island id, start and end co-ordinates, and length of the genomic island.


## Note
This program has been tested on 64-bit machine and is intended for use on 64-bit computers


