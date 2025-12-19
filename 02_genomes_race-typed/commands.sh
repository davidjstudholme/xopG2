### Download NCBI's datasets utility:
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
chmod u+x datasets

### Use NCBI's datasets utility to download the genome sequences, unzip them and make symlinks in curren directory:
./datasets download genome accession --inputfile assembly_accessions.txt --include genome --filename assemblies.zip
unzip assemblies.zip
ln -s ncbi_dataset/data/GCA_*/GCA_*.fna  ncbi_dataset/data/GCF_*/GCF_*.fna  .
rm README.md 


### Make symlinks to the genome sequence files such that symlinks have informative names and appropriate extensions for input to PhaME:
perl rename_files.pl genomes.txt



