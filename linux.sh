# nohup后台运行
nohup bash /home/shijian/project/cuttag/cuttag.sh > cuttag.log 2>&1 &

#服务器里的R
source ~/.bashrc
conda info -e
conda activate epigenetics
/home/shijian/miniconda3/envs/epigenetics/bin/R
conda activate R
/home/shijian/miniconda3/envs/R/bin/R

#查看磁盘空间大小
df -h
