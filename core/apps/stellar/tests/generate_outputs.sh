#!/bin/sh
#
# Output generation script for stellar

STELLAR=../../../../build/debug_linux/core/apps/stellar/stellar

# ============================================================
# Varying error rates
# ============================================================

eps="e-1"
errRate=0.1
${STELLAR} -d 512_simSeq1_$eps.fa -q 512_simSeq2_$eps.fa -e $errRate -l 50 -x 10 -k 7 -n 5000 -s 10000 -f -v -t -of gff -o $eps.gff > $eps.stdout

eps="75e-3"
errRate=0.075
${STELLAR} -d 512_simSeq1_$eps.fa -q 512_simSeq2_$eps.fa -e $errRate -l 50 -x 10 -k 7 -n 5000 -s 10000 -f -v -t -of gff -o $eps.gff > $eps.stdout

eps="5e-2"
errRate=0.05
${STELLAR} -d 512_simSeq1_$eps.fa -q 512_simSeq2_$eps.fa -e $errRate -l 50 -x 10 -k 7 -n 5000 -s 10000 -f -v -t -of gff -o $eps.gff > $eps.stdout

eps="25e-3"
errRate=0.025
${STELLAR} -d 512_simSeq1_$eps.fa -q 512_simSeq2_$eps.fa -e $errRate -l 50 -x 10 -k 7 -n 5000 -s 10000 -f -v -t -of gff -o $eps.gff > $eps.stdout

eps="e-4"
errRate=0.0001
${STELLAR} -d 512_simSeq1_$eps.fa -q 512_simSeq2_$eps.fa -e $errRate -l 50 -x 10 -k 7 -n 5000 -s 10000 -f -v -t -of gff -o $eps.gff > $eps.stdout


# ============================================================
# Varying minimal lengths
# ============================================================

minLen="20"
${STELLAR} -d 512_simSeq1_5e-2.fa -q 512_simSeq2_5e-2.fa -e 0.05 -l $minLen -x 10 -k 7 -n 5000 -s 10000 -f -v -t -of gff -o minLen$minLen.gff > minLen$minLen.stdout

minLen="150"
${STELLAR} -d 512_simSeq1_5e-2.fa -q 512_simSeq2_5e-2.fa -e 0.05 -l $minLen -x 10 -k 7 -n 5000 -s 10000 -f -v -t -of gff -o minLen$minLen.gff > minLen$minLen.stdout


# ============================================================
# Output format
# ============================================================

eps="5e-2"
errRate=0.05
${STELLAR} -d 512_simSeq1_$eps.fa -q 512_simSeq2_$eps.fa -e $errRate -l 50 -x 10 -k 7 -n 5000 -s 10000 -f -v -t -of text -o $eps.txt > $eps"txt.stdout"  
