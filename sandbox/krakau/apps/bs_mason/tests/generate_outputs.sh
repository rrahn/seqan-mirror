#!/bin/sh
#
# Output generation script for bs_mason

BS_MASON=../../../../../../build/S_Make/sandbox/krakau/apps/bs_mason/bs_mason

# ============================================================
# Simulate methylation rates 
# ============================================================

${BS_MASON} illumina_bs -N 100 -n 36 -s 0 -rnp read -o illumina-se-adeno-N100-n36.fasta adeno-genome.fa  > illumina-se-adeno-N100-n36.stdout
${BS_MASON} illumina_bs -N 100 -n 36 -s 0 -rnp read -o hg18_21_tiny-se-N100-n36.fasta hg18_21_tiny.fa  > hg18_21_tiny-se-N100-n36.stdout
${BS_MASON} illumina_bs -mp -N 100 -n 36 -s 0 -rnp read -o hg18_21_tiny-pe-N100-n36.fasta hg18_21_tiny.fa  > hg18_21_tiny-pe-N100-n36.stdout


# ============================================================
# Use given methylation rates
# ============================================================

${BS_MASON} illumina_bs -N 100 -n 36 -s 0 -rnp read -o hg18_21_tiny-se-N100-n36-umr.fasta -umr -mr h1_21_tiny_methLevel.txt hg18_21_tiny.fa  > hg18_21_tiny-se-N100-n36-umr.stdout
${BS_MASON} illumina_bs -mp -N 100 -n 36 -s 0 -rnp read -o hg18_21_tiny-pe-N100-n36-umr.fasta -umr -mr h1_21_tiny_methLevel.txt hg18_21_tiny.fa  > hg18_21_tiny-pe-N100-n36-umr.stdout


