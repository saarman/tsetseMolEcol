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

## The input files:
/uufs/chpc.utah.edu/common/home/saarman-group1/MullingOverStuff/glossina-full-scaffolds.vcf.gz
/uufs/chpc.utah.edu/common/home/saarman-group1/MullingOverStuff/filter_vcf.slurm

## The output files:
/uufs/chpc.utah.edu/common/home/saarman-group1/MullingOverStuff/glossina_a-scaffolds.vcf
/uufs/chpc.utah.edu/common/home/saarman-group1/MullingOverStuff/glossina_auto-scaffolds.vcf
/uufs/chpc.utah.edu/common/home/saarman-group1/MullingOverStuff/glossina_b-scaffolds.vcf
/uufs/chpc.utah.edu/common/home/saarman-group1/MullingOverStuff/glossina_c-scaffolds.vcf
/uufs/chpc.utah.edu/common/home/saarman-group1/MullingOverStuff/glossina_d-scaffolds.vcf
/uufs/chpc.utah.edu/common/home/saarman-group1/MullingOverStuff/glossina_e-scaffolds.vcf
/uufs/chpc.utah.edu/common/home/saarman-group1/MullingOverStuff/glossina_f-scaffolds.vcf
/uufs/chpc.utah.edu/common/home/saarman-group1/MullingOverStuff/glossina_na-scaffolds.vcf
/uufs/chpc.utah.edu/common/home/saarman-group1/MullingOverStuff/glossina_sex-scaffolds.vcf
