# Calculate SNV types for each variant calling algorithm;

use strict;

my ($vcf, $out, $sample)=@ARGV;
open (VCF, $vcf) or die;
open (OUT, ">>", $out) or die;

my ($AC,$AG,$AT,$CA,$CG,$CT,$GA,$GC,$GT,$TA,$TC,$TG)=(0,0,0,0,0,0,0,0,0,0,0,0);
while(<VCF>){
	chomp;
	my ($ref, $alt)=(split/\t/,$_)[3,4];
	if($ref=~/A/){
		if($alt=~/C/){$AC++;}
		if($alt=~/G/){$AG++;}
                if($alt=~/T/){$AT++;}
	}
        if($ref=~/C/){
                if($alt=~/A/){$CA++;}
                if($alt=~/G/){$CG++;}
                if($alt=~/T/){$CT++;}
        }
        if($ref=~/G/){
                if($alt=~/A/){$GA++;}
                if($alt=~/C/){$GC++;}
                if($alt=~/T/){$GT++;}
        }
        if($ref=~/T/){
                if($alt=~/A/){$TA++;}
                if($alt=~/C/){$TC++;}
                if($alt=~/G/){$TG++;}
        }
}
close(VCF);

print OUT $sample."\tA"."\tC"."\tG"."\tT"."\n";
print OUT "A"."\t"."0\t".$AC."\t".$AG."\t".$AT."\n".
	"C"."\t".$CA."\t"."0\t".$CG."\t".$CT."\n".
	"G"."\t".$GA."\t".$GC."\t"."0\t".$GT."\n".
	"T"."\t".$TA."\t".$TC."\t".$TG."\t0"."\n";
