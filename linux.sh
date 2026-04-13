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
df -h
df -h /
#清空临时文件夹
sudo rm -rf /tmp/*   cme811811
#这个是清理系统日志
sudo journalctl --vacuum-time=3d  
#清空日志文件
du -sh /var/log
rm -rf /var/log/*.log

#已删除但仍占空间
lsof | grep deleted
