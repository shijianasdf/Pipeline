#! /bin/bash
# https://hub.docker.com/ 在该网址上面查询下载已有的容器 biofly i l m
#展示镜像
docker image ls
#查找所有容器
docker ps -a 
#下载docker容器
docker pull data2intelligence/data2intelligence-suite 
#加载容器
docker run -it data2intelligence/data2intelligence-suite /bin/bash  
#加载容器并创建共享文件夹
docker run -it -v /Users/biofly/testCytosig:/data2intelligence/containerDir data2intelligence/data2intelligence-suite /bin/bash 






