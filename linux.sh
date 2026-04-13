# nohup后台运行
nohup bash /home/shijian/project/cuttag/cuttag.sh > cuttag.log 2>&1 &

#服务器里的R
source ~/.bashrc
conda info -e
conda activate epigenetics
/home/shijian/miniconda3/envs/epigenetics/bin/R
conda activate R
/home/shijian/miniconda3/envs/R/bin/R

#查看linux版本
hostnamectl
#查看磁盘结构
lsblk
#查看磁盘空间大小
df -h   #查看物理设备和磁盘挂载点，挂载到根目录下的磁盘为系统盘
df -h /
#清空临时文件夹
sudo rm -rf /tmp/*   cme811811
sudo rm -rf /var/tmp/*
#这个是清理系统日志
sudo journalctl --vacuum-time=3d  
#清空日志文件
du -sh /var/log
sudo rm -rf /var/log/*.log
#系统缓存
sudo apt clean
#已删除但仍占空间
lsof | grep deleted

# 重点排查区，排查home目录下文件大小
sudo du -h --max-depth=1 /home | sort -hr
# 找超大文件
sudo find / -type f -size +5G

#已删除但仍占空间（高级坑🔥）
sudo lsof | grep deleted

#为用户在/data目录下创造私有空间
sudo mkdir /data/xiewanhua
sudo chown xiewanhua:xiewanhua /data/xiewanhua
sudo chmod 700 /data/xiewanhua
