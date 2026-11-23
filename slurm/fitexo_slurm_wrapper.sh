/#!/bin/bash

# Получаем паттерн из аргумента командной строки
PATTERN="$1"

# Если паттерн не передан, показываем ошибку
if [ -z "$PATTERN" ]; then
    echo "Usage: $0 'pattern/*.bam'"
    echo "Example: $0 '../data/reads/D_melanogaster_mnase/SRR*/*.sorted_bam'"
    exit 1
fi

echo "Pattern: $PATTERN"

# Раскрываем паттерн в массив
files=($PATTERN)
num_files=${#files[@]}

if [ $num_files -eq 0 ]; then
    echo "Error: No files found matching: $PATTERN"
    exit 1
fi

echo "Found $num_files files:"
for i in "${!files[@]}"; do
    echo "  [$i] ${files[$i]}"
done

echo "Submitting array job 0-$((num_files-1))"

# Запускаем массив задач
sbatch --array=0-$((num_files-1)) fitexo_slurm.sh "$PATTERN"