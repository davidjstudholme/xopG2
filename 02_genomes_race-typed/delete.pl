#!/usr/bin/perl -w

use strict;
use warnings;


my $delete_list_file = 'for-deletion.txt ';
my $genomes_file = 'old-genomes.txt';

my %numbers;
open(DEL_FILE, "<$delete_list_file") or die $!;
while(<DEL_FILE>) {
    if (m/WHRI\s*(\d+)/) {
	my $number = $1;
	warn "$1\n";
	$numbers{$number}++;
    }
}


open(GEN_FILE, "<$genomes_file") or die $!;
while (<GEN_FILE>) {
    chomp;
    if (m/WHRI\s*(\d+)/) {
	if (defined $numbers{$1}) {
	    warn "Ignore line: $_\n";
	} else {
	    print "$_\n";
	}
    }
}
