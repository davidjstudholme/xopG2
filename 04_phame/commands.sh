### This takes a while to run, so best to do it in a screen session
#screen

### Already installed PhaME into a Conda environment
#conda activate phame_env
conda list -n phame_env > phame_env_packages.txt

### Assumes that PhaME is already installed
phame ./phame.fasttree.ctl

