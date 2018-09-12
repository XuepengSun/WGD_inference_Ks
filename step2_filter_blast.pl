#!/usr/bin/perl -w
use Bio::SeqIO;

my ($cdsFile,$blsatp) = @ARGV;

my $io=new Bio::SeqIO(-file=>$cdsFile,-format=>'fasta');
while($seq=$io->next_seq){
	$hash{$seq->id} = $seq->seq;
}
$io->close;

open IN,$blastp || die;
open OUT,">data_for_mcl.txt";

while(<IN>){
	my @s = split;
	if($s[0] eq $s[1]){
		print OUT $s[0],"\t",$s[1],"\t",$s[-2],"\n";
		next;
	}
	next if ($s[4] / $s[2] < 0.5 || $s[4] / $s[3] < 0.5);
	my $ck = 0;
	if($s[-2] > 98){
		open A,">seq1.fa";
		open B,">seq2.fa";
		print A ">",$s[0],"\n",$hash{$s[0]},"\n";
		print B ">",$s[1],"\n",$hash{$s[1]},"\n";
		close A;
		close B;
		system("blastn -query seq1.fa -subject seq2.fa -evalue 1e-5 -outfmt \"6 qseqid sseqid qlen slen length qstart qend sstart send pident bitscore\" -out blastn.out");
		open TMP,"blastn.out";
		while(<TMP>){
			@s=split;
			$ck++ if $s[-2] > 98;
			last;
		}
		close TMP;
	}
	if($ck == 0){
		print OUT $s[0],"\t",$s[1],"\t",$s[-2],"\n";
	}
}
close IN;
close OUT;

system("rm seq1.fa seq2.fa blastn.out");

