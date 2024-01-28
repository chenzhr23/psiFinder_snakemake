#!/usr/bin/perl
use strict;
use warnings;

if ( scalar(@ARGV) != 1 ) {
	die "Usage: perl $0 <annoFile>\n";
	die "Example: perl $0 bedAnnotator_res.bed >out 2>log\n";
}

my $annoFile  = $ARGV[0];

#### new xuan start ###
my @protein_coding = ("protein_coding");
my @pseudogene = ("rRNA_pseudogene", "transcribed_unprocessed_pseudogene", "translated_processed_pseudogene", "unprocessed_pseudogene", "processed_pseudogene", "transcribed_processed_pseudogene", "unitary_pseudogene", "pseudogene", "polymorphic_pseudogene", "transcribed_unitary_pseudogene", "TR_V_pseudogene", "TR_J_pseudogene", "IG_V_pseudogene", "IG_C_pseudogene", "IG_D_pseudogene", "IG_pseudogene", "IG_J_pseudogene");
my @lncRNA = ("lncRNA","processed_transcript","lincRNA","non_coding","3prime_overlapping_ncRNA","3prime_overlapping_ncrna","sense_intronic","antisense","sense_overlapping","known_ncrna","macro_lncRNA","bidirectional_promoter_lncRNA","retained_intron","TEC");
my @sncRNA = ("snRNA","snoRNA","misc_RNA","miscRNA","miRNA","ribozyme","sRNA","scRNA","scaRNA","srpRNA","tRNA-Deu","tRNA-RTE","piRNA","siRNA");
my @rRNA = ("rRNA","Mt_rRNA");
my @tRNA = ("tRNA","Mt_tRNA","vaultRNA");
my @IG_gene = ("IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_V_gene");
my @TR_gene = ("TR_C_gene","TR_J_gene","TR_V_gene","TR_D_gene");
my @repeatMasker = ("5S-Deu-L2","Alu","centr","CR1","DNA","DNA?","ERV1","ERV1?","ERVK","ERVL","ERVL?","ERVL-MaLR","Gypsy","Gypsy?","hAT","hAT?","hAT-Ac","hAT-Blackjack","hAT-Charlie","hAT-Tag1","hAT-Tip100","hAT-Tip100?","Helitron","Helitron?","L1","L2","Low_complexity","LTR","LTR?","MIR","MULE-MuDR","nonsense_mediated_decay","non_stop_decay","Penelope","PiggyBac","PiggyBac?","RNA","RTE-BovB","RTE-X","Satellite","Simple_repeat","SVA","TcMar?","TcMar-Mariner","TcMar-Tc2","TcMar-Tigger","telo","Unknown","acro","Crypton","Dong-R4","I-Jockey","Kolobok","L1-Tx1","Merlin","MULE-MuDR?","PIF-Harbinger","SINE?","TcMar","TcMar-Pogo","TcMar-Tc1");
my @intergenic = ("intergenic");
my @circRNA = ("circRNA");
# my $readthrough = ("readthrough", "stop_codon_readthrough", "readthrough_transcript");
#### new xuan end ###

my $outFile = $annoFile;
$outFile =~ s/\.bed/\.biotype\.bed/;
open(IN,"<$annoFile") or die "can't open $annoFile:$!\n";
open(OUT,">$outFile") or die "can't open $outFile:$!\n";
while(<IN>){
	chomp;
	if($_ !~ /^#/){
		my @lines = split("\t", $_);
		# my $gene = $lines[2];
		my $info = $lines[6];
		if($info eq "intergenic"){
			print  OUT join("\t", @lines[0..6]),"\t","intergenic","\t","intergenic","\t",join("\t", @lines[7..11]),"\n";
		}else{
			my @infoArr = split(/\|/, $info);
			my $geneType = $infoArr[4];
			# print $geneType,"\n";
			if(grep { $_ eq $geneType } @circRNA){
				print OUT join("\t", @lines[0..7]),"\t","circRNA","\t",join("\t", @lines[8..12]),"\n";
			}elsif(grep { $_ eq $geneType } @intergenic){
				print OUT join("\t", @lines[0..7]),"\t","intergenic","\t",join("\t", @lines[8..12]),"\n";
			}elsif(grep { $_ eq $geneType } @protein_coding){
				print OUT join("\t", @lines[0..7]),"\t","mRNA","\t",join("\t", @lines[8..12]),"\n";
			}elsif(grep { $_ eq $geneType } @IG_gene){
				print OUT join("\t", @lines[0..7]),"\t","mRNA","\t",join("\t", @lines[8..12]),"\n";
			}elsif(grep { $_ eq $geneType } @TR_gene){
				print OUT join("\t", @lines[0..7]),"\t","mRNA","\t",join("\t", @lines[8..12]),"\n";
			}elsif(grep { $_ eq $geneType } @lncRNA){
				print OUT join("\t", @lines[0..7]),"\t","lncRNA","\t",join("\t", @lines[8..12]),"\n";
			# }elsif(grep { $_ eq $geneType } @sncRNA){
			# 	print OUT join("\t", @lines[0..7]),"\t",$geneType,"\t",join("\t", @lines[8..12]),"\n";
			}elsif(grep { $_ eq $geneType } @sncRNA){
				print OUT join("\t", @lines[0..7]),"\t","sncRNA","\t",join("\t", @lines[8..12]),"\n";
			}elsif(grep { $_ eq $geneType } @rRNA){
				print OUT join("\t", @lines[0..7]),"\t","rRNA","\t",join("\t", @lines[8..12]),"\n";
				}elsif(grep { $_ eq $geneType } @tRNA){
				print OUT join("\t", @lines[0..7]),"\t","tRNA","\t",join("\t", @lines[8..12]),"\n";
			}elsif(grep { $_ eq $geneType } @repeatMasker){
				print OUT join("\t", @lines[0..7]),"\t","repeatMasker","\t",join("\t", @lines[8..12]),"\n";
			}elsif(grep { $_ eq $geneType } @pseudogene){
				print OUT join("\t", @lines[0..7]),"\t","pseudogene","\t",join("\t", @lines[8..12]),"\n";
			}
		}
	}
}
close(IN);
close(OUT);