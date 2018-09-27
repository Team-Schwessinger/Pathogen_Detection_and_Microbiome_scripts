# Pathogen_Detection_and_Microbiome_scripts
This repository is used for depositing the scripts and explainations supplimentary to Hu et al. 'Pathogen Detection and Microbiome Analysis of Infected Wheat Using a Portable DNA Sequencer' paper.

Here are brief descriptions for each script, detailed explainations are commented within the code:


QC_and_BLAST.py: combined all the programes used in the analysis and save ouput files from each step into the desired place.

creating_final_dataframe.py: re-organizing all the output files of basecalling and each BLAST step, summarizing all the relavant information for data visualization.

Prepare_figure2&3_info.ipynb: generating dataframe for plotting figure 2 and figure 3

figure_2&3_plotting.ipynb: used for plotting figure 2 and figure 3

Prepare_figure4_info.ipynb: generating dataframe and plotting figure 4

Prepare_figure5_info.ipynb: generating dataframe for plotting figure 5

figure5_plotting.ipynb: used for plotting figure 5

separate_stago_reads_into_fasta.ipynb: separate Parastagonospora nodorum reads from each replicates.
The Parastagonospora nodorum reads from each replicates were then concatenated together before performing BLAST analysis.

suppliment_figureS1_plotting.ipynb: ultilizing BLAST output file from the BLAST analysis of Stago reads against NCBI nt database, plotting supplimentary figureS1
