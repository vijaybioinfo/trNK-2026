#!/bin/sh
#SBATCH --job-name=R24_Cancer_Batches-1-to-20_NK_ALL_rean
#SBATCH --output=/path/to/user/sequencing_data/05-14-2022/jobs_scripts/reanalyze/R24_Cancer_Batches-1-to-20_NK_ALL/R24_Cancer_Batches-1-to-20_NK_ALL_04-01-2024_1/R24_Cancer_Batches-1-to-20_NK_ALL_cellranger_reanalyze.out.txt
#SBATCH --error=/path/to/user/sequencing_data/05-14-2022/jobs_scripts/reanalyze/R24_Cancer_Batches-1-to-20_NK_ALL/R24_Cancer_Batches-1-to-20_NK_ALL_04-01-2024_1/R24_Cancer_Batches-1-to-20_NK_ALL_cellranger_reanalyze.err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=250g
#SBATCH --time=120:00:00

echo -e '### --------------------------- Arguments --------------------------- ###'
echo -e 'Name for data analysis: R24_Cancer_Batches-1-to-20_NK_ALL'
echo -e 'General output directory: /mnt/hpcscratch/vfajardo/sequencing_data/05-14-2022/reanalyze'
echo -e 'Job output directory: /path/to/user/sequencing_data/05-14-2022/jobs_scripts/reanalyze/R24_Cancer_Batches-1-to-20_NK_ALL/R24_Cancer_Batches-1-to-20_NK_ALL_04-01-2024_1'
echo -e 'Original aggregation output: /path/to/user/sequencing_data/05-14-2022/reanalyze/data/R24_Cancer_Batches-1-to-20_NK_ALL_OrigAggr.h5'

cd /mnt/hpcscratch/vfajardo/sequencing_data/05-14-2022/reanalyze
/path/to/home/bin/cellranger-7.1.0/cellranger reanalyze --id=R24_Cancer_Batches-1-to-20_NK_ALL --matrix=/path/to/user/sequencing_data/05-14-2022/reanalyze/data/R24_Cancer_Batches-1-to-20_NK_ALL_OrigAggr.h5 --barcodes=/path/to/user/sequencing_data/05-14-2022/reanalyze/data/R24_Cancer_Batches-1-to-20_NK_ALL_BarcodeSubset.csv --localcores=2 --localmem=250 --jobmode=slurm --maxjobs=50 --mempercore=16

echo -e 'Job completed!\nCheck for errors if any.'
