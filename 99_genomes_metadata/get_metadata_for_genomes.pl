#!/usr/bin/perl

use strict ;
use warnings ;

my $genomes_list_file = './genomes.txt';

### Get a list of assembly accession numbers
my %assembly_accessions;
open my $fh, "<", $genomes_list_file or die "Cannot open file 'genomes_list_file': $!";
while (<$fh>) {
    if (m/^\#/) {
	warn "Ignore $_\n";
    } else {     
	if  (m/(GC[AF]_\d+\.\d+)/) {
	    warn $1;
	    $assembly_accessions{$1}++;
	}
    }
}

### Define which elements of assembly metadata will be queried
my @elements = qw(AssemblyAccession
		  AssemblyName
		  Organism
		  SpeciesTaxid
		  BioSampleAccn
		  ContigN50
		  Coverage);
my @elements_tmp = @elements;
my $elements_string = shift @elements_tmp;
while (my $element = shift @elements_tmp) {
    $elements_string .= ",$element";
}
warn $elements_string;

### Get a table of assemblies
my @assemblies;
foreach my $assembly_accession (sort keys %assembly_accessions) {
    warn "Assembly '$assembly_accession'\n";
    my $cmd = "esearch -db assembly -query $assembly_accession | efetch -format docsum | \
       	    xtract -pattern DocumentSummary -element $elements_string";
    warn "$cmd\n";
    open my $fh, "-|", $cmd or die "Cannot run $cmd: $!";
    while (my $readline = <$fh>) {
	chomp $readline;
	my %metadata;
	my @fields = split /\t/, $readline, -1;

	if (@fields != @elements) {
	    warn "WARNING: Field count mismatch for $assembly_accession\n";
	    warn "Expected: " . scalar(@elements) . ", got: " . scalar(@fields) . "\n";
	    warn "Line content: '$readline'\n";
	}

	foreach my $element (@elements) {
	    $metadata{$element} = shift @fields // "";
	}
	push @assemblies, \%metadata;
    }
}

### Get the strain, pathovar and owner and BioProject link associated with each assembly
foreach my $metadata_ref (@assemblies) {
    warn "Assembly: $$metadata_ref{AssemblyAccession}\n";

    ### Get strain
    my $cmd = "esearch -db biosample -query $$metadata_ref{'BioSampleAccn'} | \
       	      	   efetch -format xml | \
    		   xtract -pattern Attribute -if '\@attribute_name' -equals \"strain\" -element .";
    my $result = `$cmd`;
    chomp $result;
    #warn "\t$result\n";
    if ($result =~ m/Attribute\s+\"(.+)\"/) {
	my $strain = $1;
        warn "\tStrain=$strain\n";
	$$metadata_ref{'Strain'} = $strain;	
    }

    ### Get pathovar 
    $cmd = "esearch -db biosample -query $$metadata_ref{'BioSampleAccn'} | \
	    efetch -format xml | \ 
	    xtract -pattern Attribute -if '\@attribute_name' -equals \"pathovar\" -element .";
    $result = `$cmd`;
    chomp $result;
    #warn "\t$result\n";
    if ($result =~ m/Attribute\s+\"(.+)\"/) {
	my $pathovar = $1;
	warn "\tPathovar=$pathovar\n";
	$$metadata_ref{'Pathovar'} = $pathovar;
    }

    ### Get owner
    $cmd = "esearch -db biosample -query $$metadata_ref{'BioSampleAccn'} | \
    	   efetch -format docsum | \
    	   xtract -pattern Owner -element Name";
    $result = `$cmd`;
    chomp $result;
    my $owner = $result;
    warn "\tOwner=$owner\n";
    $$metadata_ref{'Owner'} = $owner;

    ### Get BioProject Link(s)
    $cmd = "esearch -db biosample -query $$metadata_ref{'BioSampleAccn'} | \
           efetch -format docsum | \ 
           xtract -pattern Links -element Link";
    $result = `$cmd`;
    chomp $result;
    $result =~ s/\t/ /g; # Remove any internal tab characters that might arise form there being multiple links
    my $link = $result;
    warn "\tLink=$link\n";
    $$metadata_ref{'Link'} = $link;

}


### Print header line
my @elements_tmp = (@elements, 'Strain', 'Pathovar', 'Owner', 'Link');
my $first_element = shift @elements_tmp;
print "$first_element";
foreach my $element (@elements_tmp) {
    print "\t$element";
}
print "\n";

### List the assemblies
foreach my $metadata_ref (@assemblies) {
    my @elements_tmp = (@elements, 'Strain', 'Pathovar', 'Owner', 'Link');
    my $first_element = shift @elements_tmp;
    print "$$metadata_ref{$first_element}";
    foreach my $element (@elements_tmp) {
	if (defined $$metadata_ref{$element}) {
	    print "\t$$metadata_ref{$element}";
	} else {
	    print "\t";
	}
    }
    print "\n";
}

