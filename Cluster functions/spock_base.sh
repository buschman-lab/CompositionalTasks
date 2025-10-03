
#SBATCH -o out/%x_%A_%a.out
#SBATCH -e err/%x_%A_%a.err
#SBATCH -p all
#SBATCH --mail-type FAIL
#SBATCH --mail-user tafazoli@princeton.edu


echo "In the directory: `pwd` "
echo "As the user: `whoami` "
echo "on host: `hostname` "
echo "Num Cores: ${NUM_CORES}"
echo "Array Allocation Number: $SLURM_ARRAY_JOB_ID"
echo "Array Index: $SLURM_ARRAY_TASK_ID"

module load matlab/R2021b

