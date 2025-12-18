
python3 extract_gbk.py ssrS pepQ

conda create -n clinker_env
conda activate clinker_env
conda install -c conda-forge clinker

conda list -n bwa_env > clinker_env_packages.txt
conda env export > clinker_env.yaml

clinker -p ssrS-pepQ.clinker.html  *ssrS_to_pepQ.gbk


