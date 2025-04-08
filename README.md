# tsetseMolEcol
Filter tsetse VCF file by Muller Element for Molecular Ecology Class Spring BIOL4750/6750

Make sure to allocate resources before you start
```
salloc --time=1:00:00 --ntasks 1 --mem=100G --account=saarman-np --partition=saarman-shared-np
```

## Filter by Muller Element with grep

Testing with A.txt
```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/MullingOverStuff
chmod -R g+w ./*

bash
zgrep "^#" glossina-full-scaffolds.vcf.gz > glossina-A-scaffolds.vcf.gz

# test with A.txt
while read scaffold; do
  zgrep -E "^${scaffold}\s" glossina-full-scaffolds.vcf.gz >> glossina-A-scaffolds.vcf
done < A.txt

chmod -R g+w ./*
```

Looping over each regions.txt that Eric uploaded
```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/MullingOverStuff
chmod -R g+w ./*
bash

# Grab the VCF header once
zgrep "^#" glossina-full-scaffolds.vcf.gz > vcf_header.tmp

# Loop over each scaffold list in the scaffolds/ directory
for file in scaffolds/*.txt; do
  # Extract base name for output (e.g., glossina_a.txt → glossina-a-scaffolds.vcf)
  base=$(basename "$file" .txt)
  outfile="${base}-scaffolds.vcf"

  # Start the output file with the header
  cp vcf_header.tmp "$outfile"

  # Append matching scaffold lines from the VCF
  while read scaffold; do
    zgrep -E "^${scaffold}\s" glossina-full-scaffolds.vcf.gz >> "$outfile"
  done < "$file"

  echo "Created $outfile"
done

# Clean up temp header file
rm vcf_header.tmp

bgzip *.vcf
chmod -R g+w ./*
```

## Filter with BCFtools, too many errors 
```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/MullingOverStuff
chmod -R g+w ./*

module load htslib
bgzip glossina-full-scaffolds.vcf
tabix -p vcf glossina-full-scaffolds.vcf.gz
# error = sorting is needed before tabix

bcftools sort glossina-full-scaffolds.vcf.gz -Oz -o glossina-full-scaffolds-sorted.vcf.gz
# error = no scaffold numbers in header lines

module load bcftools
bcftools view -r A.txt glossina-full-scaffolds.vcf.gz -oZ glossina-A-scaffolds.vcf.gz
# error = needs to be gz zipped
```
Where regions.txt contains one scaffold per line, or regions like:
```
scaffold_1
scaffold_2:1-3842897
```



