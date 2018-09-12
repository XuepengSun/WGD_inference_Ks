
#1 filter proteins
  perl step1_filter_protein.pl protein.fa > filtered.pep.fa

#2 self blasting
  makeblastdb -in filtered.pep.fa -dbtype prot

  ##split query file for speed
  blastp -db filtered.pep.fa -query filtered.pep.fa \
  -num_threads 64 -evalue 1e-5 \
  -outfmt "6 qseqid sseqid qlen slen length qstart qend sstart send pident bitscore" \
  -out filtered.pep.selfblast 

#3 filter blastp (output: data_for_mcl.txt)
  perl step2_filter_blast.pl cds.fasta filtered.pep.selfblast 

#4 MCL clustering
  mcl data_for_mcl.txt --abc -I 1.5 -o mcl.1.5.txt 

#5 Calculating dNdS (output: dNdS_combined_output.txt)
  perl dNdS.pl mcl.1.5.txt cds.fasta protein.fasta 

#6 plotting
  run KsPlots.R

