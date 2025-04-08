# tsetseMolEcol
Filter tsetse VCF file by Muller Element for Molecular Ecology Class Spring BIOL4750/6750

Make sure to allocate resources before you start
```
salloc --time=1:00:00 --ntasks 1 --mem=100G --account=saarman-np --partition=saarman-shared-np
```

Filter with BCFtools
```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/MullingOverStuff
chmod -R g+w ./*

module load bcftools
bcftools view -r regions.txt input vcf -o output vcf

```
Where regions.txt contains one scaffold per line, or regions like:
```
scaffold_1
scaffold_2:1-3842897
```
