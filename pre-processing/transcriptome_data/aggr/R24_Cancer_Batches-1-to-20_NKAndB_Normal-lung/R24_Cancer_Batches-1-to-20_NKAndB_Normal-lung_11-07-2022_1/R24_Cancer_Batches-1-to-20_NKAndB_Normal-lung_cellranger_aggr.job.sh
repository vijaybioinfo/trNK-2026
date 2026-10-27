#PBS -N R24_Cancer_Batches-1-to-20_NKAndB_Normal-lung_aggr
#PBS -o /mnt/hpcscratch/vfajardo/sequencing_data/05-14-2022/jobs_scripts/aggr/R24_Cancer_Batches-1-to-20_NKAndB_Normal-lung/R24_Cancer_Batches-1-to-20_NKAndB_Normal-lung_11-07-2022_1/R24_Cancer_Batches-1-to-20_NKAndB_Normal-lung_cellranger_aggr.out.txt
#PBS -e /mnt/hpcscratch/vfajardo/sequencing_data/05-14-2022/jobs_scripts/aggr/R24_Cancer_Batches-1-to-20_NKAndB_Normal-lung/R24_Cancer_Batches-1-to-20_NKAndB_Normal-lung_11-07-2022_1/R24_Cancer_Batches-1-to-20_NKAndB_Normal-lung_cellranger_aggr.err.txt
#PBS -m abe
#PBS -M vfajardo@lji.org
#PBS -q default
#PBS -l nodes=1:ppn=8
#PBS -l mem=250gb
#PBS -l walltime=80:00:00

echo -e '\n\n######### ------------- Job to Run Seurat Analysis ------------ #########\n\n'
echo -e '### --------------------------- Arguments --------------------------- ###'
echo -e '### --------------------- Job general arguments --------------------- ###'
echo -e 'Name for data aggregation: R24_Cancer_Batches-1-to-20_NKAndB_Normal-lung'
echo -e 'General output directory: /mnt/hpcscratch/vfajardo/sequencing_data/05-14-2022/aggr'
echo -e 'Job output directory: /mnt/hpcscratch/vfajardo/sequencing_data/05-14-2022/jobs_scripts/aggr/R24_Cancer_Batches-1-to-20_NKAndB_Normal-lung/R24_Cancer_Batches-1-to-20_NKAndB_Normal-lung_11-07-2022_1'
echo -e 'csv file with count matrices to be aggregated: /path/to/user/sequencing_data/05-14-2022/aggr/data/R24_Cancer_Batches-1-to-20_NKAndB_Normal-lung_aggr_table.1.0.csv'

cd /mnt/hpcscratch/vfajardo/sequencing_data/05-14-2022/aggr
/path/to/home/bin/cellranger-3.1.0/cellranger aggr --id=R24_Cancer_Batches-1-to-20_NKAndB_Normal-lung --csv=/path/to/user/sequencing_data/05-14-2022/aggr/data/R24_Cancer_Batches-1-to-20_NKAndB_Normal-lung_aggr_table.1.0.csv --normalize=mapped --localcores=8 --localmem=250 --nosecondary --jobmode=torque

echo -e 'Job completed!\nCheck for errors if any.'
