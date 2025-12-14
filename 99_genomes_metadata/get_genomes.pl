#!/usr/bin/perl

use strict ;
use warnings ;

if (@ARGV) {
    warn "Will get assemblies for these BioProjects: @ARGV\n";
} else {
    die "Usage $0 <list of BioProjects e.g. PRJNA1209959 PRJNA1040293>\n";
}

### Define the BioProjects for inclusion
my @bioprojects = @ARGV;

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
foreach my $bp (@bioprojects) {
    warn "BioProject '$bp'\n";
    my $cmd = "esearch -db assembly -query $bp | efetch -format docsum | \
       	    xtract -pattern DocumentSummary -element $elements_string";
    warn "$cmd\n";
    open my $fh, "-|", $cmd or die "Cannot run $cmd: $!";
    while (my $readline = <$fh>) {
	chomp $readline;
	my %metadata;
	$metadata{'BioProject'} = $bp;
	my @fields = split /\t/, $readline, -1;

	if (@fields != @elements) {
	    warn "WARNING: Field count mismatch for BioProject $bp\n";
	    warn "Expected: " . scalar(@elements) . ", got: " . scalar(@fields) . "\n";
	    warn "Line content: '$readline'\n";
	}

	foreach my $element (@elements) {
	    $metadata{$element} = shift @fields // "";
	}
	push @assemblies, \%metadata;
    }
}

### Get the strain, pathovar and owner  associated with each assembly
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
}


### Print header line
my @elements_tmp = ('BioProject', @elements, 'Strain', 'Pathovar', 'Owner');
my $first_element = shift @elements_tmp;
print "$first_element";
foreach my $element (@elements_tmp) {
    print "\t$element";
}
print "\n";

### List the assemblies
foreach my $metadata_ref (@assemblies) {
    my @elements_tmp = ('BioProject', @elements, 'Strain', 'Pathovar', 'Owner');
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

