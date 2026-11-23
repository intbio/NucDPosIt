#!/bin/bash

#SBATCH --job-name fitexo
#SBATCH --partition intbio
#SBATCH --cpus-per-task 8
#SBATCH --mem 8G
#SBATCH --output slurm_out_%A_%a.log
#SBATCH --error ./slurm_errors_%A_%a.log
#SBATCH --time=01:00:00

# Получаем паттерн из аргумента
PATTERN="$1"

# Раскрываем паттерн в массив
files=($PATTERN)

# Берем файл по индексу задачи
file_path="${files[$SLURM_ARRAY_TASK_ID]}"

echo "========================================="
echo "Array Job ID: $SLURM_ARRAY_JOB_ID"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURMD_NODENAME"
echo "File: $file_path"
echo "Start time: $(date)"
echo "========================================="

# Запуск обработки
conda run -n newNucDPosIT python3 $FITEXOPATH \
    -s 0.0000005 -e 0.00030 -st 0.00000025 \
    -n 3 -@ $SLURM_CPUS_PER_TASK -i "$file_path" -f

echo "Finished at: $(date)"
