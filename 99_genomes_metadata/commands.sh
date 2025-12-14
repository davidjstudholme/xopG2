### Install NCBI EDirect tool
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
export PATH=${HOME}/edirect:${PATH}

./get_genomes.pl PRJNA1209959 PRJNA1040293 > test.tsv

./get_genomes.pl PRJNA1209959 PRJNA1040293 PRJNA991952 PRJNA689092 PRJNA185979 PRJNA185977 PRJNA991952 PRJNA991733 PRJNA729255 PRJNA689092 PRJNA641237 PRJNA185980  PRJNA63187  PRJNA1040293 PRJNA185979 PRJNA185977 PRJNA774128 PRJNA988133 PRJNA505826 PRJNA774128 PRJNA485808 > genomes.tsv

./get_metadata_for_genomes.pl > genomes-metadata.tsv


### Illustrative commands:

### Generate a table of assemblies
esearch -db assembly -query PRJNA991733 | efetch -format docsum | \
    xtract -pattern DocumentSummary -element AssemblyAccession,AssemblyName,Organism,SpeciesTaxid,WGS,BioSampleAccn,ContigN50,Coverage

### Get strain name for BioSamples
esearch -db biosample -query SAMN36315931 | efetch -format xml | \
    xtract -pattern Attribute -if '@attribute_name' -equals "strain" -element .

### Get pathovar for BioSamples
esearch -db biosample -query SAMN36315931 | efetch -format xml | \
    xtract -pattern Attribute -if '@attribute_name' -equals "pathovar" -element .

### Get owner of BioSample
esearch -db biosample -query SAMN36315931 | efetch -format docsum | \
    xtract -pattern Owner -element Name



