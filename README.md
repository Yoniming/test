# ACO  

  ACO_ align 用于对齐的.fastq文件的压缩，ACO_unalign用于非对齐的.fastq文件的压缩。  
  
  在这两个文件夹下，都有两个名称相同的子文件夹，encoder用于压缩，decoder用于解压。
## 使用指南
在linux操作系统下，运行程序仅需以下三个步骤：  
  1.克隆  
  
  执行命令 `git clone https://github.com/Yoniming/ACO.git` 将源码克隆到当前路径。  
  
  2.编译  
  
  在`ACO_ align`和`ACO_unalign`的`./encoder`目录下分别执行`g++ -o aco.out *.cpp` 生成可执行文件`aco.out`  
  
  在`ACO_ align`和`ACO_unalign`的`./decoder`目录下分别执行`g++ -o de_aco.out *.cpp` 生成可执行文件`de_aco.out` 
  
  3.运行
  
  将××.fastq文件放在与`aco.out`同级目录下：  
  * 若是对齐的.fastq文件，放在`../ACO_ align/encoder`下，使用命令`aco.out <文件名> [压缩模式]` 进行压缩，如:  
  `aco.out ERR233152.fastq 2`  
  * 若是非对齐的.fastq文件，放在`../ACO_ unalign/encoder`下，使用命令`aco.out <文件名>` 进行压缩，如:  
  `aco.out ERR233152.fastq`  
  
  运行完后，会在当前目录下生成一个`.stream`后缀的压缩文件，若要解压该文件，将`.stream`和`.fasta`放在对应的`./decoder` 文件夹下：  
  使用命令`de_aco.out <文件名>` 进行解压，如：  
  `de_aco.out ERR233152`

