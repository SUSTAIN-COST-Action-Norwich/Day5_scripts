#!/usr/bin/perl -w
#
#
use warnings;
use strict;
use Bio::SeqIO;
#
#
#
#--------------------------------------------------------------------------
# bioperl_make_secretome.pl
#--------------------------------------------------------------------------
#
# 
# 
#
# Usage: perl bioperl_make_secretome.pl <fasta_input_file> <targetp_output_file> <tmhmm_output_file> <fasta_output_file> <output_tab_seperated_results_file>
#
#
# Uses BioPerl
#
#
# This script parses outputfiles from "targetp" and "tmhmm" programs and filter 
# signalp2 predicted secreted proteins for presence of transmembrane domain and
# mitochondrial targeting signal. If any present, the sequence is omitted from
# the final secretome prediction.
#
# @joewinnz
# TSL 2015
#
#
#--------------------------------------------------------------------------
#

my $usage= 'Usage: perl bioperl_make_secretome.pl <fasta_input_file> <targetp_output_file> <tmhmm_output_file> <fasta_output_file> <output_tab_seperated_results_file>';
my ($fasta_infile, $targetp_file, $tmhmm_file, $fasta_outfile, $out_results) = @ARGV;

unless ($fasta_infile && $targetp_file && $tmhmm_file && $fasta_outfile && $out_results) {
	die "$usage\nAll parameters must be supplied!\n";
}
unless (-e $fasta_infile) {die "$fasta_infile doesn\'t exist\!\n"}
unless (-e $tmhmm_file) {die "$tmhmm_file doesn\'t exist\!\n"}
unless (-e $targetp_file) {die "$targetp_file doesn\'t exist\!\n"}


# hash to store all info from "tmhmm" and "targetp" outputs
my %annotations = ();

# processes the output from "tmhmm" program
open (TMHMM, $tmhmm_file) || die "Cannot open $tmhmm_file: $!\n";
while (my $line = <TMHMM>) {
	chomp $line;
	if ($line =~ /^\s+$/){
		next;
	}
	my @columns = split ("\t", $line);
	my $id = $columns[0];
	my $helix_info = $columns[4];
	my $helices = 0;
	if ($helix_info =~ /^PredHel=(\d+)/) {
		$helices = $1;
	} else {
		print "The output file $tmhmm_file is not from TMHMM2.0c\n";
		exit;
	}
	my ($len) = $columns[1] =~ /len=(\S+)/;
	my ($expAA) = $columns[2] =~ /ExpAA=(\S+)/;
	my ($first60) = $columns[3] =~ /First60=(\S+)/;
	my ($predHel) = $columns[4] =~ /PredHel=(\S+)/;
	my ($topology) = $columns[5] =~ /Topology=(\S+)/;
	$annotations{$id}{'helices'} = $helices;
	$annotations{$id}{'tmhmm'} = "$len\t$expAA\t$first60\t$predHel\t$topology";
	#print "$id has $helices helices\n";
}
close TMHMM;

# processes the output from "targetp" program
open (TARGETP, $targetp_file) || die "Cannot open $targetp_file: $!\n";
while (<TARGETP>) {
	chomp;
	next if (/^$/);
	next if (/^\s+$/);
	next if (/^-/);
	my @columns = split (" +");
	next unless (($columns[1] =~ /\d+/) && ($columns[1] > 0));
	my $id = $columns[0];
	my $localization = $columns[5];
	#print "$id\t$localization\n";
	$annotations{$id}{'localization'} = $localization;
	my $data = join("\t",splice @columns,2,7);
	$annotations{$id}{'targetp'}=$data;
}
close TARGETP;


open (OUT,">$out_results") || die "Cannot create output file: $out_results: $!";
print OUT "ID\tSP hmm score\tCleavage Site\tProtein length\tExpAA\tFirst60\tPredHel\tTopology\tmTP\tSP\tOther\tLoc\tRC\tTP length\tSequence\n";

my $input_fasta = Bio::SeqIO->new (-file => $fasta_infile,
							-format => 'fasta');

my $output_fasta = Bio::SeqIO->new (	-file =>">$fasta_outfile",
								-format => 'fasta');
my $mito_count = 0;	
my $pex_count = 0;
my $seq_count = 0;
my $tm_count = 0;
my $chloro_count = 0;

# processes the input sequences and collect relevant information from annotation hash
while (my $seq_obj = $input_fasta->next_seq) {
	$seq_count++;
	my $seq_id = $seq_obj->display_id;
	if (exists $annotations{$seq_id}) {
		if ($annotations{$seq_id}{'helices'} > 0) {
		# discard the sequence if it contains one or more trans-membrane domain/helix
			$tm_count++;
			if ($annotations{$seq_id}{'localization'} eq "M") {
				$mito_count++;
			}
			if ($annotations{$seq_id}{'localization'} eq "C") {
				$chloro_count++;
			}
			next;
		}
		if (($annotations{$seq_id}{'localization'} eq "S")||($annotations{$seq_id}{'localization'} eq "_")) {
			my $descr = $seq_obj->desc;
			my ($sp_hmm) = $descr =~ /--HMM=(\S+)/;
			my ($cleavage_site) = $descr =~ /--NNCleavageSite=(\S+)/;
			$seq_obj->desc($descr." --TargetP=$annotations{$seq_id}{'localization'} --TM=0");
			$output_fasta->write_seq($seq_obj);
			my $seq = $seq_obj->seq;
			print OUT "$seq_id\t$sp_hmm\t$cleavage_site\t$annotations{$seq_id}{'tmhmm'}\t$annotations{$seq_id}{'targetp'}\t$seq\n";
			$pex_count++;
		} else {
			# discard the sequence if it is predicted to contain either mitochondrial
			# or chloroplastic targeting signal
			if ($annotations{$seq_id}{'localization'} eq "M") {
				$mito_count++;
			}
			if ($annotations{$seq_id}{'localization'} eq "C") {
				$chloro_count++;
			}
		}
	}
}
close OUT;

print "\nTotal number of sequences processed: $seq_count\n";
print "Total number of sequences with at least one TM: $tm_count\n";
print "Total number of sequences with SP, no TM, and predicted secreted: $pex_count\n";
print "Total number of sequences with mitochondrial signal: $mito_count\n";
print "Total number of sequences with chloroplast signal: $chloro_count\n";



exit;

