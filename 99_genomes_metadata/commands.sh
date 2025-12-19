### Install NCBI EDirect tool
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
export PATH=${HOME}/edirect:${PATH}

./get_metadata_for_genomes.pl > genomes-metadata.tsv





