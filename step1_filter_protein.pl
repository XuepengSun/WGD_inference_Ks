#!/usr/bin/perl -w
use Bio::SeqIO;

my $maxProteinSize = 3000;
my $minProteinSize = 100;

my $io = new Bio::SeqIO(-file=>$ARGV[0],-format=>'fasta');

while($seq=$io->next_seq){
	my $xl = $seq->seq;
	$xl=~s/\*$//;
	# remove premature sequence
	next if $xl =~/\*/;
	next if $seq->length < $minProteinSize;
	next if $seq->length > $maxProteinSize;
	print ">",$seq->id,"\n",$seq->seq,"\n";
}
$io->close;


