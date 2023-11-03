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


 * In cases of competitive exclusion (species i always drives species j to extinction), we inferred that 𝛼𝑖𝑗 < 1 and 𝛼𝑗𝑖 > 1. For bistability (the high-abundance species drives the low-abundance one to extinction), we inferred that 𝛼𝑖𝑗 > 1 and 𝛼𝑗𝑖 > １．

Commensialism – where one species benefits while the other is unaffected.
Mutualism – both species benefit.
Parasitism – one species benefits while one is harmed.
Competition – neither benefits.
Predation – one species benefits while the other dies, and
Neutralism – both species unaffected.
Amensalism is an ecological interaction between two species, but in this association among organisms of two different species, one is destroyed or inhibited, and the other remains unaffected

||<1|=1|>1|
|<1|Mutualism|Commensialism|Predation|
|=1|Commensialism|Neutralism|Amensalism|
|>1|Predation/Parasitism|amensalism|Competition|

| |<1|＝１｜>１｜
｜－｜－|－｜－--df｜
｜<１｜ｃｏｏｐｅｒａｔｉｏｎ｜｜｜
｜＝１｜｜｜｜
｜>１｜｜｜ｃｏｍｐｅｔｉｔｉｏｎ｜



Numerical methods
We modeled the long-term dynamics and diversity of ecological communities using the wellknown generalized Lotka-Volterra (gLV) model, modified to include dispersal from a species pool:

