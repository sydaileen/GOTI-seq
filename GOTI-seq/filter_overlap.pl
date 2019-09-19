# filter variants not overlapped with variants called from other methods;

use strict;

my ($in1, $in2, $in3, $out)=@ARGV;
open (IN1, $in1) or die;
open (IN2, $in2) or die;
open (IN3, $in3) or die;
open (OUT, ">", $out) or die;

my %var1;
foreach(<IN1>){
	chomp;
	if(/^#/){
		next;
	}else{
		my @array=(split/\t/,$_);
		$var1{$array[0]."\t".$array[1]}="1";
	}
}
close(IN1);

my %var2;
foreach(<IN2>){
        chomp;
        if(/^#/){
                next;
        }else{
                my @array=(split/\t/,$_);
                $var2{$array[0]."\t".$array[1]}="1";
        }
}
close(IN2);


while(<IN3>){
	chomp;
	if(/^#/){
		next;
	}else{
		my $line=$_;
		my ($chr,$pos)=(split/\t/,$_)[0,1];
		if(exists $var1{$chr."\t".$pos}){
			if(exists $var2{$chr."\t".$pos}){
				print OUT $line."\n";
			}
		}
	}
}
close(IN3);
