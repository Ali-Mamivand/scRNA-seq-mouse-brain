# scRNA-seq mouse brain
In this study, I have analyzed an example of Single-cell RNA sequencing data using combination of different R packages.
The datasets analyzed during current study are available in GEO(GSE118068) which belongs to the article “Childhood cerebellar 
tumors mirror conserved fetal transcriptional programs”. Due to the limitations in my computational recourses, only the two 
timepoints E10 and E12 have been analyzed.
Please follow the next steps: 
1.	First of all, download GSM3317999(E10) and GSM3318000(E12) datasets from GEO.
	You should download three files for each dataset including genes, barcodes and matrix (genes.tsv.gz, 
	barcodes.tsv.gz, and matrix.mtx.gz). Please make sure you have downloaded all the files. 
2.	Create a new directory and put the code in this folder (Working directory).
3.	Create two new folders in the working directory called E10 and E12.
4.	Run scRNAseq.R in the working directory.
