#下载安装miniconda
wget https://mirrors.ustc.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh
wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-Latest-Linux-x86_64.sh
wget --no-check-certificate https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
mv Miniconda3-latest-Linux-x86_64.sh /pub6/temp/shijian/
cd /pub6/temp/shijian
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc　

#来查看已经安装的软件
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
conda install sra-tools=2.9.6 -y
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
conda install -c hcc aspera-cli #conda安装aspera
conda install meme homer
conda install kallisto
conda install -c bioconda cutadapt


#删除miniconda
rm -rf /home/shijian2015/miniconda3

#下载安装aspera
wget https://download.asperasoft.com/download/sw/connect/3.9.1/ibm-aspera-connect-3.9.1.171801-linux-g2.12-64.tar.gz
tar zxvf ibm-aspera-connect-3.9.1.171801-linux-g2.12-64.tar.gz
bash ibm-aspera-connect-3.9.1.171801-linux-g2.12-64.sh

#加入环境路径
echo 'export PATH=$PATH:~/.aspera/connect/bin' >> ~/.bash_profile
source ~/.bash_profile
