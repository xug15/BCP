# BCP
Bacteria Community Predictor

## Process steps.

* Abundence data to calculate interaction
* Aplication sequence data to annotation pathway.
* Using machine learning to train model.

## Annotation gene and pathway

### 1. download software and install.
• 数据库下载地址如下： 
https://github.com/microDM/MicFunPred

### 2.Set abundance matrix
```r
rowname=paste("Species",1:48,sep='')
colname=paste("s",1:48,sep='')
dyn48=data.frame(diag(48))
colnames(dyn48)=colname
rownames(dyn48)=rowname
write.table(dyn48,'d:/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/16s.abun.tsv',sep="\t",quote=F)
```
注意,文件中,16s.abun.tsv第一行表头需要加一个制表符.
 

### 3. run annotation
```sh
ln -s /mnt/d/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/16s_sequence.fa otu.fa
ln -s /mnt/d/a-document/sequencing_center_desktop/yanglab/dynamic_project/science.abm7841_data_s1/16s.abun.tsv otu.abun.txt

abun='otu.abun.txt'
fasta='otu.fa'
output='bacteria_function'
echo "MicFunPred_run_pipeline.py -i ${abun} -r ${fasta} -o ${output} -t 10 --verbose"
MicFunPred_run_pipeline.py -i ${abun} -r ${fasta} -o ${output} -t 10 --verbose
```


