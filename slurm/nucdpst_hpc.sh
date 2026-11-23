#!/bin/bash
#SBATCH --job-name=yeasts_genes
#SBATCH --output=logs/em_model_%A_%a.out
#SBATCH --error=logs/em_model_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=intbio
#SBATCH --gres=gpu:nvidia_rtx_a5000:1                
#SBATCH --time=24:0:0


conda run -n newNucDPosIT \
    python3 ~/NucDPosIT/scripts/nucdposit.py \
    ~/NucDPosIT/data/reads/S_cerevisiae/SRR1802185/SRR1802185.sorted_bam \
    ~/NucDPosIT/data/reads/S_cerevisiae/SRR1802185/exomodel_results/errors.csv \
    --regions_path "$HOME/NucDPosIT/data/s_cerevisiae_genes/s_cerevisiae_genes_1000extended.bed"\
    --output "$HOME/NucDPosIT/data/genome_scan/SRR1802185_genes.csv" \
    --device cuda \
    --min_iter 100 \
    --max_iter 1000 \
    --dyad_dist 100 \
    --nfits 100