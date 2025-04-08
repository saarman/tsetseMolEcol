#!/bin/bash

#SBATCH --partition=saarman-shared-np             # Partition to run the job
#SBATCH --account=saarman-np                      # Account to charge job resources
#SBATCH --time=24:00:00                           # Maximum runtime (24 hours)
#SBATCH --mem=24576                               # Memory in MB
#SBATCH --nodes=1                                 # Number of nodes
#SBATCH --ntasks-per-node=4                       # Number of CPU cores per node
#SBATCH --job-name="zgrep"                          # Job name for SLURM queue

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
