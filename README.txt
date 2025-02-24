Parallel implementation of set cover problem apply to gene fusion 


Requirements:
	Linux x64
	MPI
	A processing cluster
	A SAM data file



In order to launch the program, follow the steps bellow:

1 - Organize files
Copy the program files wherever you want. Be careful, you need to have the permission to edit some files. Having the administration privilegies may help.
Create your output directory, and all the necessary other directories.

2 - Place the script folder
The script directory can stay with the source code, it's not really necessary to move it but you can if necessary.

3 - Download the refGene
The version of the file we use is in the script directory but you can download/create an another one if you want.
You can go to https://genome.ucsc.edu/cgi-bin/hgTables in order to generate one.
If you do this, you need to have the same columns in the same order than our file (just download it without filter columns).

4 - Config
You need to configure some files:
	- code/makefile
		You may have to change MPI's version or some dependency but be careful on what you modify, this could lead the program to a non working state.
	- init_ref_genes_files.sh
		The script parameters depend on how you organize your files
	- script/chrNames.txt
		In this file:
			the left column corresponds to the denomination of the chromosomes in the refGene file you downloaded
			the right column corresponds to the denomination of the chromosomes in your SAM header
	- GeneFusion.sh
		In this file, every variable in the Parameter section should be redefined.
		The variables names are explicit.

5 - GeneFusion
	Basic use:
		Connect to the cluster
		time bash GeneFusion.sh
	Advance use options:
		If you pass one parameter(n) to the script, it won't ask you if you made the configuration.
		If you pass a second parameter(y) to the script, it will automatically remove the existing files in your output directory.
