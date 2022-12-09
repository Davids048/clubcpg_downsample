# clubcpg_downsample
Using downsample to account for errors caused by read depth differences among samples

NOTE: every step is performed on CluBCpG output on one chromosome

## 1. Create cluster patterns
* The cluster patterns in the clubcpg_cluster output are not consistent among different bins. So before performing downsampling, the user need to use this script to create cluster pattern labels that are consistent.
* This script takes a list of pathnames (absolute) of clubcpg outputs, the output directory, and the file name (start with "/" and ends with ".csv")
	* The list of path names of clubcpg outputs can be of any length, so that the user can put as many clubcpg-cluster output as they want to make the labels consistent across multiple analysis.
	
Sample command:
```
 python get_cluster_pattern.py --samples "clubcpg output directory 1" "clubcpg output directory 2" -o "output directory" --name "/filename.csv"
```

## 2. Create least common read depth
* Create a .csv file summarizing the lowest read depth seen in each bin
* The list of path names of clubcpg outputs can be of any length, so that the user can put as many clubcpg-cluster output as they want to make the least common read depth consistent across multiple analysis.

	Sample command:
	```
 	python get_lc_reads.py --samples "clubcpg output directory 1" "clubcpg output directory 2" -o "output directory" --name "/filename.csv"
	```

## 3. Perform downsampling procedure
* use the previously generated two files and clubcpg output to perform downsampling.
You can run ``` python downsample.py --help``` to see the commands
	* required input:
		* Cluster pattern path
		* lowest common read depth path
		* -A : clubcpg output path
		* If a user want to perform downsampling on two samples from 2 different clubcpg outputs (e.g.: downsample on sample A from clubcpg output 1 and sample B from clubcpg output 2), they have the following optional arguments:
			* -B : path of clubcpg output 2
			* -sampleA: identify which sample to choose from clubcpg output 1, allowed inputs are "A" or "B"
			* -sampleB: identify which sample to choose from clubcpg output 2, allowed inputs are "A" or "B"
			* NOTE: to perform downsampling on two samples from 2 different clubcpg outputs, the user must use the three arguments above.
	* recommended input:
		* -ncore: the number of cpu cores that will be used for running this script. Using more cores will speed up the process.
		* -o : the output directory
		* -name: the name of the output file (start with "\" and ends with .csv)
		Sample command:
		```
		 python downsample.py "cluster_patterns path" "lowest common read depths path" -A "clubcpg output path" -B "clubcpg output path 2" -ncore 10 -o "output directory path" -name "/output2.csv"
		```


