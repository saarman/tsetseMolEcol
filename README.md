# tsetseMolEcol
Filter tsetse VCF file by Muller Element for Molecular Ecology Class Spring BIOL4750/6750


Filter with BCFtools
```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/MullingOverStuff
chmod -R g+w ./*

module load bcftools
bcftools view -r A.txt input vcf -o output vcf

```
