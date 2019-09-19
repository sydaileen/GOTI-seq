# awk variation info from annotated file of ANNOVAR

use strict;
my ($in, $out)=@ARGV;
open (IN, $in) or die;
open (OUT, ">", $out) or die;

while (<IN>){
	chomp;
	if (/^Chr/){
		next;
	}
	my @line=(split/\t/,$_);
	my $type="c.".$line[1].$line[3].">".$line[4];
	my @mut=(split/:/,$line[-2])[3,4,-2,-1];
	my @wt=(split/:/,$line[-1])[3,4,-2,-1];
	my $mut_alt=$mut[0]+$mut[1];
	my $mut_ref=$mut[2]+$mut[3];
	my $wt_alt=$wt[0]+$wt[1];
	my $wt_ref=$wt[2]+$wt[3];
	my $altfreq=$mut_alt/($mut_alt+$mut_ref);
	if($line[5]=~/exonic/){
		my $pro=(split/:/,$line[9])[-1];
		$type=$pro."/".$type;
	}
	print OUT $line[0]."\t".$line[1]."\t".$type."\t".$line[5]."\t".$line[6]."\t".$mut_alt."\t".$mut_ref."\t".$wt_alt."\t".$wt_ref."\t".$altfreq."\n";
}
close(IN);
