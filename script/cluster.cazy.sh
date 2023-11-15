#!/bin/bash
#SBATCH -o /public/home/2022122/xugang/condaenv/BCP/script/cazy_1.%j.out
#SBATCH -e /public/home/2022122/xugang/condaenv/BCP/script/cazy_1.%j.error
#SBATCH --partition=Gnode
#SBATCH -J carzy
#SBATCH -N 1
#SBATCH -n 28

source /public/home/2022122/xugang/bashrc

conda run -n py38 python learn_code_meaning.py