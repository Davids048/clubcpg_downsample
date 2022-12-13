# clubcpg_downsample
Using downsample to account for errors caused by read depth differences among samples
NOTE: every step is performed on CluBCpG output on one chromosome
NOTE: all pathname used in the scripts are absolute path
* Downsample procedure flowchart
	![Overview drawio](https://user-images.githubusercontent.com/90978028/207406162-a3da2561-58b7-4704-a37a-8e92b736cb5b.png)
* The procedure uses output from `clubcpg-cluster` command in the package CluBCpG. For more information about CluBCpG, you can read the documentation [here](https://clubcpg.readthedocs.io/en/latest/index.html)
* The three blue blocks represent the three scripts that are used in the procedure. And the yellow block represent the procedure of downsampling and the requierd files.


## 1. Create cluster patterns
* The cluster patterns in the clubcpg_cluster output are not consistent among different bins. So before performing downsampling, the user need to use this script to create cluster pattern labels that are consistent.
* Inputs:
	* `--samples`: a list of pathnames to clubcpg-cluster output files (.csv), 
	* `-o`: the output directory
	* `--name`: the file name (start with "/" and ends with ".csv")
	* NOTE: The list of path names of clubcpg outputs can be of any length, so that the user can put as many clubcpg-cluster output as they want to make the labels consistent across multiple analysis.
* Output: a .csv file containing the labels to the methylation patterns found in the given list of clubcpg-cluster outputs.
	
Sample command:
```
 python get_cluster_pattern.py --samples "clubcpg output directory 1" "clubcpg output directory 2" -o "output directory" --name "/filename.csv"
```

## 2. Create lowest common read depth
* Create a .csv file summarizing the lowest read depth seen in each bin
* Inputs:
	* `--samples`: a list of pathnames to clubcpg-cluster output files (.csv), 
	* `-o`: the output directory
	* `--name`: the file name (start with "/" and ends with ".csv")
	* NOTE: The list of path names of clubcpg outputs can be of any length, so that the user can put as many clubcpg-cluster output as they want to make the labels consistent across multiple analysis.
* Output: a .csv file containing the lowest commmon read depth of each bin in the given list of clubcpg-cluster outputs.

	Sample command:
	```
 	python get_lc_reads.py --samples "clubcpg output directory 1" "clubcpg output directory 2" -o "output directory" --name "/filename.csv"
	```

## 3. Perform downsampling procedure
* use the previously generated two files and clubcpg output to perform downsampling.
You can run ``` python downsample.py --help``` to see the commands
	* required input:
		* pathname to the created cluster pattern file (.csv file)
		* pathname to the created lowest common read depth file (.csv file).
		* -A : pathname to a selected clubcpg-cluster output that contains the sample of interest for downsampling.
		* If a user want to perform downsampling on two samples from 2 different clubcpg outputs (e.g.: downsample on sample A from clubcpg output 1 and sample B from clubcpg output 2), they have the following optional arguments:
			* -B :pathname to the 2nd clubcpg-cluster output file that contains the 2nd sample of interest for downsampling.
			* -sampleA: identify which sample to choose from clubcpg output 1 (-A), allowed inputs are "A" or "B"
			* -sampleB: identify which sample to choose from clubcpg output 2 (-B), allowed inputs are "A" or "B"
			* NOTE: to perform downsampling on two samples from 2 different clubcpg outputs, the user must use the three arguments above.
			* For more information on how to use 2 file mode, see section 4.
	* recommended input:
		* -ncore: the number of cpu cores that will be used for running this script. Using more cores will speed up the process.
		* -o : the output directory
		* -name: the name of the output file (start with "\" and ends with .csv)
		Sample command:
		```
		 python downsample.py "cluster_patterns path" "lowest common read depths path" -A "clubcpg output path" -B "clubcpg output path 2" -ncore 10 -o "output directory path" -name "/output2.csv"
		```

## 4. Example for using downsampling:
* **Background**: Assume we want to compare **male, neuron, p12 (12 days after birth)** data and **male,neuron,p35** data.
	* If what we have is clubcpg-cluster output files from an analysis of the two desired samples (i.e. a cluster output file comparing reads from **male, neuron, p12 (12 days after birth)** data and **male,neuron,p35** data
		* Then we only need one sample, so the arguments needed is only -A
	* If the desired samples are in two separate clubcpg-cluster output data, then we need to use 2 file mode;
		* Assume we have 2 cluster files: 
			1. clubcpg-cluster output 1 on **male, neuron, p12** and **female neuron p12**, where in the "class split" column, A ~male and B~female
			2. clubcpg-cluster output 2 on **male, neuron, p35** and **female neuron p35**, where in the "class split" column, A ~male and B~female
		* Then we need the following arguments:
			* -A : pathname to clubcpg-cluster output 1
			* -B : pathname to clubcpg-cluster output 2
			* --sampleA  "A"
			* --sampleB  "B"
* Summary of example
* 	<img width="819" alt="image" src="https://user-images.githubusercontent.com/90978028/207411373-8db5c3ad-5734-448b-902e-18f377ef77e3.png">


