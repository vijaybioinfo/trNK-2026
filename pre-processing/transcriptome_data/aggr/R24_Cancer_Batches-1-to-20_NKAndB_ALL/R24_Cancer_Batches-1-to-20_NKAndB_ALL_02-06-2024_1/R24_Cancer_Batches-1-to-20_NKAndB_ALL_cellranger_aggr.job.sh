#!/bin/sh
#SBATCH --job-name=R24_Cancer_Batches-1-to-20_NKAndB_ALL_aggr
#SBATCH --output=/mnt/hpcscratch/vfajardo/sequencing_data/05-14-2022/jobs_scripts/aggr/R24_Cancer_Batches-1-to-20_NKAndB_ALL/R24_Cancer_Batches-1-to-20_NKAndB_ALL_02-06-2024_1/R24_Cancer_Batches-1-to-20_NKAndB_ALL_cellranger_aggr.out.txt
#SBATCH --error=/mnt/hpcscratch/vfajardo/sequencing_data/05-14-2022/jobs_scripts/aggr/R24_Cancer_Batches-1-to-20_NKAndB_ALL/R24_Cancer_Batches-1-to-20_NKAndB_ALL_02-06-2024_1/R24_Cancer_Batches-1-to-20_NKAndB_ALL_cellranger_aggr.err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250g
#SBATCH --time=60:00:00

echo -e '\n\n######### ------------- Job to Run Seurat Analysis ------------ #########\n\n'
echo -e '### --------------------------- Arguments --------------------------- ###'
echo -e '### --------------------- Job general arguments --------------------- ###'
echo -e 'Name for data aggregation: R24_Cancer_Batches-1-to-20_NKAndB_ALL'
echo -e 'General output directory: /mnt/hpcscratch/vfajardo/sequencing_data/05-14-2022/aggr'
echo -e 'Job output directory: /mnt/hpcscratch/vfajardo/sequencing_data/05-14-2022/jobs_scripts/aggr/R24_Cancer_Batches-1-to-20_NKAndB_ALL/R24_Cancer_Batches-1-to-20_NKAndB_ALL_02-06-2024_1'
echo -e 'csv file with count matrices to be aggregated: /path/to/user/sequencing_data/05-14-2022/aggr/data/R24_Cancer_Batches-1-to-20_NKAndB_ALL_aggr_table.1.0.csv'

cd /mnt/hpcscratch/vfajardo/sequencing_data/05-14-2022/aggr
/path/to/home/bin/cellranger-3.1.0/cellranger aggr --id=R24_Cancer_Batches-1-to-20_NKAndB_ALL --csv=/path/to/user/sequencing_data/05-14-2022/aggr/data/R24_Cancer_Batches-1-to-20_NKAndB_ALL_aggr_table.1.0.csv --normalize=mapped --localcores=8 --localmem=250 --nosecondary --jobmode=slurm --maxjobs=20 --mempercore=16

echo -e 'Job completed!\nCheck for errors if any.'
