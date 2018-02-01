# CAFE

## Input
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

## Output

