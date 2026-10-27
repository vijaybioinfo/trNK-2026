#PBS -N R24E_Hu_6B_CD8_4D_CITE_HTODem
#PBS -o /path/to/user/R24/jobs_scripts/deconvolution/HTO_based/05-22-2021/R24_E/R24E_Hu_6B_CD8_4D_CITE/R24E_Hu_6B_CD8_4D_CITE_06-11-2021_100umis_thold/R24E_Hu_6B_CD8_4D_CITE_HTO_demultiplexing.out.txt
#PBS -e /path/to/user/R24/jobs_scripts/deconvolution/HTO_based/05-22-2021/R24_E/R24E_Hu_6B_CD8_4D_CITE/R24E_Hu_6B_CD8_4D_CITE_06-11-2021_100umis_thold/R24E_Hu_6B_CD8_4D_CITE_HTO_demultiplexing.err.txt
#PBS -m abe
#PBS -M vfajardo@lji.org
#PBS -q default
#PBS -l nodes=1:ppn=2
#PBS -l mem=30gb
#PBS -l walltime=02:00:00

echo -e "\n\n######### --------- Job to Run HTO-based deconvolution --------- #########\n\n"
echo -e "### --------------------------- Arguments --------------------------- ###"
echo -e "### --------------------- Job general arguments --------------------- ###"
echo -e "Project name: R24E_Hu_6B_CD8_4D_CITE"
echo -e "Job output directory: /path/to/user/R24/jobs_scripts/deconvolution/HTO_based/05-22-2021/R24_E/R24E_Hu_6B_CD8_4D_CITE/R24E_Hu_6B_CD8_4D_CITE_06-11-2021_100umis_thold"
echo -e "### --------------- HTO-based deconvolution arguments --------------- ###"
echo -e "Absolute path to reports directory: /path/to/user/R24/deconvolution/HTO_based/05-22-2021/R24_E/R24E_Hu_6B_CD8_4D_CITE/R24E_Hu_6B_CD8_4D_CITE_06-11-2021_100umis_thold"
echo -e "Absolute path to HTO data: /path/to/user/sequencing_data/05-22-2021/count/R24_E/R24E_Hu_6B_CD8_4D_CITE/outs/raw_feature_bc_matrix"
echo -e "Gene expression data or NULL value: "
echo -e "UMI threshold for raw data filtering or NULL value: 100"

Rscript /path/to/home/scripts/deconvolution/HTO_based/HTO_based_demultiplexing.2.4.R --PrjName R24E_Hu_6B_CD8_4D_CITE --ReportsPath /path/to/user/R24/deconvolution/HTO_based/05-22-2021/R24_E/R24E_Hu_6B_CD8_4D_CITE/R24E_Hu_6B_CD8_4D_CITE_06-11-2021_100umis_thold --HTOData /path/to/user/sequencing_data/05-22-2021/count/R24_E/R24E_Hu_6B_CD8_4D_CITE/outs/raw_feature_bc_matrix --HTONames NULL --FCThold 3  --UMIThold 100 --DoReclass doublet

echo -e "Job completed!\nCheck for errors if any."
