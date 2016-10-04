#! /bin/bash -l
#SBATCH -A b2013064
#SBATCH -n 16 -p node
#SBATCH -t 10:00:00
#SBATCH -J eb_database
#SBATCH -e /proj/b2013064/private/ACTIVE/erik/GA_16_14.10X_Chromium_evaluation/LOGS/stderr.createDB.txt
#SBATCH -o /proj/b2013064/private/ACTIVE/erik/GA_16_14.10X_Chromium_evaluation/LOGS/stdout.createDB.txt
#SBATCH --mail-type=All
#SBATCH --mail-user=erik.borgstrom@scilifelab.se

echo "$(date) Running on: $(hostname)"

cd /proj/b2013064/private/ACTIVE/erik/GA_16_14.10X_Chromium_evaluation/
workon py2.7
python initial.py

echo "$(date) AllDone"
echo "$(date) AllDone" >&2
