# *** Using species sensitivity distributions to assess the toxicity of pesticides on nitrifiers and broader microbial groups: from single species tests to soil amplicon sequencing ***

### By Maria Kolovou <sup>1,2+</sup>, Eleftheria Bachtsevani <sup>3+</sup>, Fotios Bekris <sup>2</sup>, Alexandre Pedrinho <sup>2</sup>, Graeme W. Nicol <sup>3</sup>, Christina Hazard <sup>3</sup>, Evangelia S. Papadopoulou <sup>1*</sup>, Dimitrios G. Karpouzas <sup>2</sup>

### (\* corr. author)
### (\+ contributed equally to this work)

<sup>1</sup> Laboratory of Environmental Microbiology and Virology, Department of Environmental Sciences, University of Thessaly, Larissa, Greece

<sup>2</sup> Laboratory of Plant and Environmental Biotechnology, Department of Biochemistry and Biotechnology, University of Thessalyy, Larissa, Greece

<sup>3</sup> Université Claude Bernard Lyon 1, CNRS, INRAE, VetAgro Sup, Laboratoire d'Ecologie Microbienne, Villeurbanne, 69622, France 

## The provided material includes the code used in the statistical analysis of the study.

For obtaining the code the users need to open a terminal and having the [GitHub tools](https://github.com/git-guides/install-git), git-clone or download the repository, and enter the base folder. E.g:

```
$ git clone https://github.com/mariaklv04/Metsulfuron-methyl_toxicity_on_soil_microbiome.git
```

In the case of the computational methods, with the "Metsulfuron-methyl_toxicity_on_soil_microbiome" folder as working directory, and assuming that the necessary software and R packages are installed, the used code can be executed as described in this Readme.md file. The necessary datasets for performing all sequencing based analysis can be downloaded implementing the code provided in the corresponding repository folders as explained below.

## Description of the order of executed scripts.

For Fungi, Bacteria, Ammonia oxidizing bacteria and Ammonia oxidizing archea files, steps 0-2 concern the data retrieval from NCBI and preprocessing (demultiplex) and phyloseq object construction, while step 3 and the subfolders concern the actual data analysis.

0) First, it is necessary to download the sequencing data.
To do so, you need to enter the "0.DownloadData" subfolder of "Fungi", "Bacteria", "AOA" and "AOB" folders accordingly and execute the "fetch_data.sh" bash script for batch (01-02), this assumes that you are located at the working directory "Metsulfuron-methyl_toxicity_on_soil_microbiome". The NCBI submitted amplicons are includes at those batch/files.The script is based on the SRR accession numbers for each batch file and can be found in the 0.DownloadData folder as a.txt file.
Once the download is done, you need to combine all forward reads to a single file and all reverse reads to another file as well.
```
for i in {01-02}
do
	cd Fungi/0.DownloadData/batch${i}
	sh -x fetch_data.sh
	cat *_1.fastq | gzip > forward.fastq.gz
	cat *_2.fastq | gzip > reverse.fastq.gz
	cd ../../../
	cd Bacteria/0.DownloadData/batch${i}
	sh -x fetch_data.sh
	cat *_1.fastq | gzip > forward.fastq.gz
	cat *_2.fastq | gzip > reverse.fastq.gz
	cd ../../../
    cd AOA/0.DownloadData/batch${i}
	sh -x fetch_data.sh
	cat *_1.fastq | gzip > forward.fastq.gz
	cat *_2.fastq | gzip > reverse.fastq.gz
	cd ../../../
    cd AOB/0.DownloadData/batch${i}
	sh -x fetch_data.sh
	cat *_1.fastq | gzip > forward.fastq.gz
	cat *_2.fastq | gzip > reverse.fastq.gz
	cd ../../../
done
```

1) Then you need to demultiplex the data according to our own demultiplexing method using our in-house script.
This requires Flexbar v3.0.3 to be installed, and the mapping file (map_file) accordingly for Fungi, Bacteria, AOA and AOB that is also provided in each file.
A detailed description of our in-house multiplexing approach is provided in our [previous work] (https://github.com/SotiriosVasileiadis/mconsort_tbz_degr#16s).
You need to enter the folder Fungi (or Bacteria, AOA, AOB)/1.Demultiplex and run the following commands (change the MY_PROCS variable to whatever number of logical processors you have available and want to devote),
the following commands are going to save the demultiplexed files in the Fungi(or Bacteria, AOA, AOB)/1.Demultiplex/demux_out folder.
```
MY_WORKING_DIR_BASE=`pwd`
for i in {01-02}
do
  cd Fungi/1.Demultiplex
  MY_PROCS=56
  bash DemuxOwnBCsys_absPATH.sh demux_out${i} ${MY_WORKING_DIR_BASE}/Fungi/0.DownloadData/batch${i}/forward.fastq.gz ${MY_WORKING_DIR_BASE}/Fungi/0.DownloadData/batch${i}/reverse.fastq.gz fun${i}Fungi_map_file.txt ${MY_PROCS}
  cd demux_out${i}/analysis_ready
  gunzip *.gz # unzips files skipped by the Demux script
  cd ../../../../
  cd Bacteria/1.Demultiplex
  MY_PROCS=56
  bash DemuxOwnBCsys_absPATH.sh demux_ou${i} ${MY_WORKING_DIR_BASE}/Fungi/0.DownloadData/batch${i}/forward.fastq.gz ${MY_WORKING_DIR_BASE}/Fungi/0.DownloadData/batch${i}/reverse.fastq.gz bac${i}Bacteria_map_file.txt ${MY_PROCS}
  cd demux_out${i}/analysis_ready
  gunzip *.gz # unzips files skipped by the Demux script
  cd ../../../../
 cd AOA/1.Demultiplex
  MY_PROCS=56
  bash DemuxOwnBCsys_absPATH.sh demux_out${i} ${MY_WORKING_DIR_BASE}/Fungi/0.DownloadData/batch${i}/forward.fastq.gz ${MY_WORKING_DIR_BASE}/AOA/0.DownloadData/batch${i}/reverse.fastq.gz fun${i}AOA_map_file.txt ${MY_PROCS}
  cd demux_out${i}/analysis_ready
  gunzip *.gz # unzips files skipped by the Demux script
  cd ../../../../
cd AOB/1.Demultiplex
  MY_PROCS=56
  bash DemuxOwnBCsys_absPATH.sh demux_out${i} ${MY_WORKING_DIR_BASE}/Fungi/0.DownloadData/batch${i}/forward.fastq.gz ${MY_WORKING_DIR_BASE}/AOA/0.DownloadData/batch${i}/reverse.fastq.gz fun${i}AOB_map_file.txt ${MY_PROCS}
  cd demux_out${i}/analysis_ready
  gunzip *.gz # unzips files skipped by the Demux script
  cd ../../../../
done

cd Fungi/1.Demultiplex
mkdir -p demux_out/analysis_ready
cp demux_out[0-9]/analysis_ready/*.fastq demux_out/analysis_ready/
cd ../../

cd Bacteria/1.Demultiplex
mkdir -p demux_out/analysis_ready
cp demux_out[0-9]/analysis_ready/*.fastq demux_out/analysis_ready/
cd ../../
cd AOA/1.Demultiplex
mkdir -p demux_out/analysis_ready
cp demux_out[0-9]/analysis_ready/*.fastq demux_out/analysis_ready/
cd ../../
cd AOB/1.Demultiplex
mkdir -p demux_out/analysis_ready
cp demux_out[0-9]/analysis_ready/*.fastq demux_out/analysis_ready/
cd ../../
```
2) Following, the "Quality-Classification-Phyloseq Object.R" script of the Fungi(or Bacteria, AOA, AOB)/2.PhyloseqObjectPreparation folder is run in order to prepare the final phyloseq object to be used in the data analysis described below. Before running the script make sure that the necessary reference databases are found in the same folder. The taxonomic annotations of the resulting fungal and bacterial ASVs were performed using the UNITE ITS v.8.2 (04.02.2020) (Morrison-Whittle et al., 2017) and the Silva v.138 (Yilmaz et al., 2014) databases as references respectively. AOB and AOA amoA amplicons were compared with the databases of Abell et al. (2012) and Alves et al. (2018), respectively, to classify the generated ASVs taxonomically. The sample data informations is also needed for the final construction of phyloseq objects which is also included in the files accordingly as samdf.txt. 
```
cd Fungi/2.PhyloseqObjectPreparation
# fetch the databases
wget https://files.plutof.ut.ee/public/orig/1D/B9/1DB95C8AC0A80108BECAF1162D761A8D379AF43E2A4295A3EF353DD1632B645B.gz
# run the R script
Fungi Quality-Classification-Phyloseq Object.r
cd ../../
cd Bacteria/2.PhyloseqObjectPreparation
# fetch the databases
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz
tar vxf *.gz
# run the R script
Bacteria Quality-Classification-Phyloseq Object.r
cd ../../
cd AOA/2.PhyloseqObjectPreparation
# fetch the databases
wget https://files.plutof.ut.ee/public/orig/1D/B9/1DB95C8AC0A80108BECAF1162D761A8D379AF43E2A4295A3EF353DD1632B645B.gz
# run the R script
AOA Quality-Classification-Phyloseq Object.r
cd ../../
cd AOB/2.PhyloseqObjectPreparation
# fetch the databases
wget https://files.plutof.ut.ee/public/orig/1D/B9/1DB95C8AC0A80108BECAF1162D761A8D379AF43E2A4295A3EF353DD1632B645B.gz
# run the R script
AOB Quality-Classification-Phyloseq Object.r
cd ../../
```
3) Data analysis folder include subfolders for each analysis graphs supplied at the researched article "Using species sensitivity distributions to assess the toxicity of pesticides on nitrifiers and broader microbial groups: from single species tests to soil amplicon sequencing". Subfolders contain the R script to be executed for "Fungi", "Bacteria", "AOA" and "AOB" for both main and supplementary figures. Follow the scripts to replicate the results.
```

3a. Run the NMDS Analysis (Main Figure)
3b. Run the PERMANOVA and the Pair-Wise PERMANOVA Analysis (Main Figure)
3c. Run the SSDs Analysis (Main Figure)


Further on continue with the supplementary figures

3d. Run the α-diversity nanalysis (Supplementary)
3e. Run the NMDS Analysis (Supplementary)
3f. Run the PERMANOVA and the Pair-Wise PERMANOVA Analysis (Supplementary)
3g. Run the SSDs Analysis (Supplementary)

```



## Code Usage disclaimer<a name="disclaimer"></a>

The following is the disclaimer that applies to all scripts, functions, one-liners, etc. This disclaimer supersedes any disclaimer included in any script, function, one-liner, etc.

You running this script/function means you will not blame the author(s) if this breaks your stuff. This script/function is provided **AS IS** without warranty of any kind. Author(s) disclaim all implied warranties including, without limitation, any implied warranties of merchantability or of fitness for a particular purpose. The entire risk arising out of the use or performance of the sample scripts and documentation remains with you. In no event shall author(s) be held liable for any damages whatsoever (including, without limitation, damages for loss of business profits, business interruption, loss of business information, or other pecuniary loss) arising out of the use of or inability to use the script or documentation. Neither this script/function, nor any part of it other than those parts that are explicitly copied from others, may be republished without author(s) express written permission. Author(s) retain the right to alter this disclaimer at any time. This disclaimer was copied from a version of the disclaimer published by other authors in https://ucunleashed.com/code-disclaimer and may be amended as needed in the future.
