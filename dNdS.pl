#!/usr/bin/perl -w
use File::Find;
use Bio::SeqIO;
use threads;

my $maxGroup = 4;  # maximum members in the orthogroup;
my ($orthMCL,$nucl,$prot) = @ARGV;

#system("source ~/.bashrc");

my (%hash,%ck,%data);
my $io = new Bio::SeqIO(-file=>$nucl,-format=>'fasta');
while($seq=$io->next_seq){$hash{'nucl'}{$seq->id}=$seq->seq}
$io->close;

$io = new Bio::SeqIO(-file=>$prot,-format=>'fasta');
while($seq=$io->next_seq){$hash{'prot'}{$seq->id}=$seq->seq}
$io->close;

unless(-d 'PAML'){mkdir 'PAML'}

my $n = 0;
open IN,$orthMCL || die;
while(<IN>){
	chomp;
	my @s = split;
	next if (scalar (@s) == 1 || scalar (@s) > $maxGroup);
	for($i=0;$i<=$#s-1;$i++){
		for($j=$i+1;$j<=$#s;$j++){
			$n++;
			$data{$n} = [$s[$i],$s[$j]];
		}
	}
}
close IN;

print STDERR "\ntotal comparison: ",scalar keys %hash,"\n";

foreach $d(sort {$a<=>$b} keys %data){
	&fromOrthMCL($d,$data{$d}->[0],$data{$d}->[1]);
}

open OUT,">dNdS_combined_output.txt";
print OUT "ID\tNG86(dn|ds|w)\tYN00(dn|ds|w)\tLPB93(dn|ds|w)\tML_F3x4(dn|ds|w)\n";

my @ynf = glob("PAML/*.yn00");
my @mlf = glob("PAML/*.codeml");

foreach(@ynf){
	my ($name) = $_=~/PAML\/(.*)\.yn00/;
	open IN,$_ || die;
	while(<IN>){
		if(/Nei\s+\&\s+Gojobori\s+1986\.\s+dN\/dS\s+\(dN\,\s+dS\)/){
			for(1..5){$ff = <IN>}
			($w,$dn,$ds)=$ff=~/^\S+\s+(\S+)\s+\((\S+)\s+(\S+)\)/;
			$res{$name}{'NG'} = [$dn,$ds,$w];
		}
		elsif(/\(B\)\s+Yang\s+\&\s+Nielsen\s+\(2000\)\s+method/){
			for(1..8){$ff = <IN>}
			$ff=~s/^\s+|\s+$//g;
			my @t = split(/\s+/,$ff);
			$res{$name}{'YN'} = [$t[7],$t[10],$t[6]];
		}
		elsif(/^LPB93\:\s+/){
			($ds,$dn,$w)=$_=~/dS\s+=(.*)dN\s+=(.*)w\s+=(.*)/;
			$ds=~s/^\s+|\s+$//g;
			$dn=~s/^\s+|\s+$//g;
			$w=~s/^\s+|\s+$//g;
			$res{$name}{'LPB'} = [$dn,$ds,$w];
		}
		else{next}	
	}
	close IN;
}

foreach(@mlf){
	my ($name) = $_=~/PAML\/(.*)\.codeml/;
	open IN,$_ || die;
	while(<IN>){
		if(/pairwise\s+comparison\,\s+codon\s+frequencies\:\s+F3x4/){
			for(1..7){$ff = <IN>}
			($w,$dn,$ds)=$ff=~/dN\/dS=(.*)dN\s+=(.*)dS\s+=(.*)/;
                        $ds=~s/^\s+|\s+$//g;
                        $dn=~s/^\s+|\s+$//g;
                        $w=~s/^\s+|\s+$//g;
                        $res{$name}{'ML'} = [$dn,$ds,$w];
		}
	}
	close IN;
}

foreach(sort {$a<=>$b} keys %res){
	print OUT $_,"\t",$res{$_}{'NG'}->[0],'|',$res{$_}{'NG'}->[1],'|',$res{$_}{'NG'}->[2],"\t";
	print OUT $res{$_}{'YN'}->[0],'|',$res{$_}{'YN'}->[1],'|',$res{$_}{'YN'}->[2],"\t";
	print OUT $res{$_}{'LPB'}->[0],'|',$res{$_}{'LPB'}->[1],'|',$res{$_}{'LPB'}->[2],"\t";
	print OUT $res{$_}{'ML'}->[0],'|',$res{$_}{'ML'}->[1],'|',$res{$_}{'ML'}->[2],"\n";
}
close OUT;


sub fromOrthMCL{
		my ($m,$s1,$s2) = @_;
		my $prefix = "PAML/".$m;
		my $seqNucl = $prefix.'.nucl';
		my $seqProt = $prefix.'.prot';
		my $aln = $prefix.'.mafft';
		my $paml = $prefix.'.paml';
		my $yn = $prefix.'.yn00';		
		open OUT,">$seqNucl";
		print OUT '>','seq1',"\n",$hash{'nucl'}{$s1},"\n";
		print OUT '>','seq2',"\n",$hash{'nucl'}{$s2},"\n";
		close OUT;
		open OUT,">$seqProt";
		print OUT '>','seq1',"\n",$hash{'prot'}{$s1},"\n";
		print OUT '>','seq2',"\n",$hash{'prot'}{$s2},"\n";
	        close OUT;
		system("mafft --maxiterate 1000 --localpair --quiet $seqProt > $aln");
		system("perl ~/tools/paml4.9h/pal2nal.v14/pal2nal.pl $aln $seqNucl -output paml > $paml");
		open IN,'/home/sxp/tools/paml4.9h/yn00.ctl';
		open YN,">yn00.ctl";
		while(<IN>){
			$_=~s/seqfile\s+=\s+\S+/seqfile = $paml/;
			$_=~s/outfile\s+=\s+\S+/outfile = $yn/;
			print YN $_;
		}
		close IN;
		close YN;
		system("yn00 yn00.ctl");
		system("rm 2YN.* rst rst1 rub");
		
		my $codeml = $prefix.'.codeml';
		open IN,"/home/sxp/tools/paml4.9h/codeml.ctl";
		open CD,">codeml.ctl";
		while(<IN>){
			$_=~s/seqfile\s+=\s+\S+/seqfile = $paml/;
			$_=~s/outfile\s+=\s+\S+/outfile = $codeml/;
			$_=~s/treefile.*?\*/treefile = \*/;
			$_=~s/noisy\s+=.*?\*/noisy = 0 \*/;
			$_=~s/verbose\s+=.*?\*/verbose = 0 \*/;
			$_=~s/runmode\s+=.*?\*/runmode = -2 \*/;
			$_=~s/seqtype\s+=.*?\*/seqtype = 1 \*/;
			$_=~s/CodonFreq\s+=.*?\*/CodonFreq = 2 \*/;
			$_=~s/ndata\s+=\s+\d+/ndata = 1 /;
			$_=~s/model\s+=\s+\d+/model = 0 /;
			$_=~s/NSsites\s+=\s+\d+/NSsites = 0 /;
			print CD $_;
		}
		close IN;
		close CD;
		system("codeml codeml.ctl");
		system("rm 2ML.* 2NG.* rst rst1 rub");		
		system("rm $seqNucl $seqProt $aln $paml");
}

