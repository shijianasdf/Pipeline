#下载安装miniconda
wget https://mirrors.ustc.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh
wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-Latest-Linux-x86_64.sh
wget --no-check-certificate https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
mv Miniconda3-latest-Linux-x86_64.sh /pub6/temp/shijian/
cd /pub6/temp/shijian
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc　
#下载安装3.10版本miniconda（稳定版）
wget https://repo.anaconda.com/miniconda/Miniconda3-py310_24.1.2-0-Linux-x86_64.sh
bash Miniconda3-py310_24.1.2-0-Linux-x86_64.sh
source ~/.bashrc

#来查看已经安装的软件
conda --version
conda list
#添加channels
conda config --add channels conda-forge　　　
conda config --add channels defaults　　　
conda config --add channels r　　　
conda config --add channels bioconda
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/msys2/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
conda config --add channels https://mirrors.ustc.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/bioconda/
conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/bioconda/
conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/main/
conda config --set show_channel_urls yes
#删除channels,vim打开.condark文件
vim ~/.condark
#查看已经添加的channels
conda config --get channels
#移除所有channals
conda config --remove-key channels
#创建名为bioinfo的环境
conda create -y --name bioinfo python=3 
conda create -y --name bioinfo python=2
#激活bioinfo环境
source activate bioinfo
#查看系统中已有的环境
conda info -e 
conda env list
#退出环境
conda deactivate
#删除某个环境(名为bioinfo的环境)
conda remove -n bioinfo --all
#更新miniconda
conda update conda
#安装软件
conda search 软件名
conda install 软件名=版本号
conda update 软件名
conda remove 软件名

conda install sra-tools
conda install sra-tools=3.2.1 -y
conda install r-base=4.5.1 -y
conda install parallel-fastq-dump
conda install fastqc
conda install multiqc
conda install trim-galore
conda install bwa
conda install bowtie2
conda install star
conda install hisat2
conda install rsem
conda install salmon
conda install htseq
conda install bedtools deeptools
conda install bedops macs2 samtools
conda install gatk gatk4 picard
conda install subread
conda install strelka 
conda install -c hcc aspera-cli #conda安装aspera-cli，不含密匙，需要下载Aspera Connect来获取密匙
conda install meme homer
conda install kallisto
conda install -c bioconda cutadapt

#服务器用conda配置R环境（各种版本）
#R 4.3.3 是目前最稳妥、最推荐的版本，用来做单细胞、ChIP-seq、RNA-seq 都没有问题。
#如果你想“尝鲜”，可以装 R 4.4.x，但要注意可能会遇到依赖兼容性问题。
conda search r-base #首先查看都有哪些版本的R
conda create -p /data/shijian/software/envs/R_env_4.5.1 r-base=4.5.1 #用conda创建并安装指定版本和环境位置的R
conda activate /data/shijian/software/envs/R_env_4.5.1 #激活该conda环境
conda create -p /data/shijian/software/envs/R_env_4.3.3 r-base=4.3.3
conda activate /data/shijian/software/envs/R_env_4.3.3
conda install -c conda-forge r-essentials #安装R包
conda install -c conda-forge r-tidyverse 
R

#删除miniconda
rm -rf /home/shijian/miniconda3
rm -rf ~/.conda ~/.condarc ~/.continuum
nano ~/.bashrc  #找到并删除如下自动添加的初始化段落（通常靠近底部）：


#下载安装aspera-connect
wget https://download.asperasoft.com/download/sw/connect/3.9.1/ibm-aspera-connect-3.9.1.171801-linux-g2.12-64.tar.gz
tar zxvf ibm-aspera-connect-3.9.1.171801-linux-g2.12-64.tar.gz
bash ibm-aspera-connect-3.9.1.171801-linux-g2.12-64.sh
#查找ascp和ascp密匙
find ~ -type f -name "ascp"
find ~ -type f -name "asperaweb_id_dsa.openssh"
#加入环境路径
echo 'export PATH=$PATH:~/.aspera/connect/bin' >> ~/.bash_profile
source ~/.bash_profile

#安装cell ranger
tar -xzvf ~/software/cellranger-arc-2.0.2.tar.gz
tar -xzvf ~/software/cellranger-atac-2.2.0.tar.gz
tar -xzvf ~/software/cellranger-9.0.1.tar.gz
tar -xzvf ~/refData/refdata-gex-GRCm39-2024-A.tar.gz
tar -xzvf ~/refData/refdata-gex-GRCh38-2024-A.tar.gz
tar -xzvf ~/refData/refdata-cellranger-arc-GRCm39-2024-A.tar.gz
tar -xzvf ~/refData/refdata-cellranger-arc-GRCh38-2024-A.tar.gz
export PATH=/home/shijian/software/cellranger-9.0.1:$PATH
export PATH=/home/shijian/software/cellranger-atac-2.2.0:$PATH
export PATH=/home/shijian/software/cellranger-arc-2.0.2:$PATH
