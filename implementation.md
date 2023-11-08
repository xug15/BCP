# Implementation

## 1. Calculate the interaction between species. 
```r
Rscript script/a1.calculate_alpamatrix.with.plot.r data/abandence/b2.hn.b123.bacteria.abundance.csv  data/interaction/b2.hn.b123
Rscript script/a1.calculate_alpamatrix.with.plot.r data/abandence/b2.hn.b101112.bacteria.abundance.csv data/interaction/b2.hn.b101112
Rscript script/a1.calculate_alpamatrix.with.plot.r data/abandence/b2.hn.b131415.bacteria.abundance.csv data/interaction/b2.hn.b131415
Rscript script/a1.calculate_alpamatrix.with.plot.r data/abandence/b2.hn.b161718.bacteria.abundance.csv data/interaction/b2.hn.b161718
Rscript script/a1.calculate_alpamatrix.with.plot.r data/abandence/b2.hn.b192021.bacteria.abundance.csv data/interaction/b2.hn.b192021
Rscript script/a1.calculate_alpamatrix.with.plot.r data/abandence/b2.hn.b222324.bacteria.abundance.csv data/interaction/b2.hn.b222324
Rscript script/a1.calculate_alpamatrix.with.plot.r data/abandence/b2.hn.b456.bacteria.abundance.csv data/interaction/b2.hn.b456
Rscript script/a1.calculate_alpamatrix.with.plot.r data/abandence/b2.hn.b789.bacteria.abundance.csv data/interaction/b2.hn.b789
Rscript script/a1.calculate_alpamatrix.with.plot.r data/abandence/b1.ln.b101112.bacteria.abundance.csv data/interaction/b1.ln.b101112
Rscript script/a1.calculate_alpamatrix.with.plot.r data/abandence/b1.ln.b123.bacteria.abundance.csv data/interaction/b1.ln.b123
Rscript script/a1.calculate_alpamatrix.with.plot.r data/abandence/b1.ln.b131415.bacteria.abundance.csv data/interaction/b1.ln.b1314115
Rscript script/a1.calculate_alpamatrix.with.plot.r data/abandence/b1.ln.b161718.bacteria.abundance.csv data/interaction/b1.ln.b161718
Rscript script/a1.calculate_alpamatrix.with.plot.r data/abandence/b1.ln.b192021.bacteria.abundance.csv data/interaction/b1.ln.b192021
Rscript script/a1.calculate_alpamatrix.with.plot.r data/abandence/b1.ln.b222324.bacteria.abundance.csv data/interaction/b1.ln.b222324
Rscript script/a1.calculate_alpamatrix.with.plot.r data/abandence/b1.ln.b456.bacteria.abundance.csv data/interaction/b1.ln.b456
Rscript script/a1.calculate_alpamatrix.with.plot.r data/abandence/b1.ln.b789.bacteria.abundance.csv data/interaction/b1.ln.b789

```
## 2. merge interaction table.
```sh
cd data/interaction/
cat b1*interaction.csv > b3.ln.interaction.csv
grep -v '""' b3.ln.interaction.csv > b3.ln.csv
cut -f 2,3,4 -d ',' b3.ln.csv > b3.ln.interaction.csv
rm b3.ln.csv
cat b2*interaction.csv > b4.hn.interaction.csv
grep -v '""' b4.hn.interaction.csv > b4.hn.csv
cut -f 2,3,4 -d ',' b4.hn.csv > b4.hn.interaction.csv
rm b4.hn.csv
cd ../..
```

## 3. Annotation pathway and gene in bacteria.
```sh

abun='data/annotation/16s.abun.tsv'
fasta='data/annotation/16s_sequence.fa'
output='data/annotation/mic'
echo "MicFunPred_run_pipeline.py -i ${abun} -r ${fasta} -o ${output} -t 10 --verbose"
MicFunPred_run_pipeline.py -i ${abun} -r ${fasta} -o ${output} -t 10 --verbose

cp ${output}/

```
## 4. extract useful message from annotation files.
```sh
cp data/annotation/mic/out.blast data/annotation/species.txt
cp data/annotation/mic/CAZymes_metagenome/CAZymes_metagenome.tsv.gz data/annotation
cp data/annotation/mic/COG_metagenome/COG_metagenome.tsv.gz data/annotation
cp data/annotation/mic/KO_metagenome/KEGG_pathways_MinPath_prunned.tsv.gz data/annotation
cp data/annotation/mic/KO_metagenome/KO_metagenome_MinPath_prunned.tsv.gz  data/annotation
cp data/annotation/mic/MetaCyc_metagenome/Pathway_summarize_by_Types.tsv.gz  data/annotation/cyc.tsv.gz
cp data/annotation/mic/Pfam_metagenome/Pfam_metagenome.tsv.gz  data/annotation
cp data/annotation/mic/TIGRFAM_metagenome/TIGRFAM_metagenome.tsv.gz  data/annotation

gunzip  data/annotation/CAZymes_metagenome.tsv.gz
gunzip data/annotation/COG_metagenome.tsv.gz
gunzip data/annotation/KEGG_pathways_MinPath_prunned.tsv.gz
gunzip data/annotation/KO_metagenome_MinPath_prunned.tsv.gz
gunzip data/annotation/cyc.tsv.gz
gunzip data/annotation/Pfam_metagenome.tsv.gz
gunzip data/annotation/TIGRFAM_metagenome.tsv.gz

```
## 5. get species name. 
```sh
cut -f 1,2  data/annotation/species.txt > data/annotation/species_name.txt
```

## 6. combine pathway with interaction.
```sh
Rscript script/a2.combine.gene.interaction.r data/annotation/CAZymes_metagenome.tsv data/interaction/b3.ln.interaction.csv data/ml/ln.ml.cazymes 
Rscript script/a2.combine.gene.interaction.r data/annotation/COG_metagenome.tsv data/interaction/b3.ln.interaction.csv data/ml/ln.ml.cog 
Rscript script/a2.combine.gene.interaction.r data/annotation/KEGG_pathways_MinPath_prunned.tsv data/interaction/b3.ln.interaction.csv data/ml/ln.ml.kegg
Rscript script/a2.combine.gene.interaction.r data/annotation/KO_metagenome_MinPath_prunned.tsv data/interaction/b3.ln.interaction.csv data/ml/ln.ml.ko 
Rscript script/a2.combine.gene.interaction.r data/annotation/Pfam_metagenome.tsv data/interaction/b3.ln.interaction.csv data/ml/ln.ml.pfam 
Rscript script/a2.combine.gene.interaction.r data/annotation/TIGRFAM_metagenome.tsv data/interaction/b3.ln.interaction.csv data/ml/ln.ml.tigrfam 
#
Rscript script/a2.combine.gene.interaction.r data/annotation/CAZymes_metagenome.tsv data/interaction/b4.hn.interaction.csv data/ml/hn.ml.cazymes 
Rscript script/a2.combine.gene.interaction.r data/annotation/COG_metagenome.tsv data/interaction/b4.hn.interaction.csv data/ml/hn.ml.cog 
Rscript script/a2.combine.gene.interaction.r data/annotation/KEGG_pathways_MinPath_prunned.tsv data/interaction/b4.hn.interaction.csv data/ml/hn.ml.kegg
Rscript script/a2.combine.gene.interaction.r data/annotation/KO_metagenome_MinPath_prunned.tsv data/interaction/b4.hn.interaction.csv data/ml/hn.ml.ko 
Rscript script/a2.combine.gene.interaction.r data/annotation/Pfam_metagenome.tsv data/interaction/b4.hn.interaction.csv data/ml/hn.ml.pfam 
Rscript script/a2.combine.gene.interaction.r data/annotation/TIGRFAM_metagenome.tsv data/interaction/b4.hn.interaction.csv data/ml/hn.ml.tigrfam 
```



