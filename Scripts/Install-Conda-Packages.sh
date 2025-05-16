

# 1. This is folder in which we downloaded the conda 
mkdir miniconda


# 2. Download the conda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# 3. Refreshing the terminal 
source ~/miniconda3/bin/activate

# You should see the terminal that has (base)
# Like this one "(base) rafat@XB0-P17-P006-02:~/miniconda$ "

# 4. Create the conda environment in which the packages we shall use 
conda create -n RNA-seq

# 5. Activate the Environment 
conda activate RNA-seq



# 6. Installing the packages needed 
conda install -c bioconda kallisto fastqc multiqc sra-too
ls cutadapt
