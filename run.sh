abun='data/otu/ph_abun.tsv'
fasta='data/otu/pH_preferences_ASVs_PAN_ID.607.fasta'
output='data/otu/ph_mic'
#python script/otu2abun.py -i data/otu/16s_sequence.fa -o data/otu/16s_abuns.tsv  
echo "MicFunPred_run_pipeline.py -i ${abun} -r ${fasta} -o ${output} -t 10 --verbose"
MicFunPred_run_pipeline.py -i ${abun} -r ${fasta} -o ${output} -t 10 --verbose







