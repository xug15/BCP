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
output='data/annotation/bacteria_function'
echo "MicFunPred_run_pipeline.py -i ${abun} -r ${fasta} -o ${output} -t 10 --verbose"
MicFunPred_run_pipeline.py -i ${abun} -r ${fasta} -o ${output} -t 10 --verbose


```




