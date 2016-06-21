Genome Coverage Suite
===
This program allows users to graph coverage from chip sequencing. There are five main features, getting coverage data from bam files, subtracting a control, graphing the whole genome, graphing a range of positions in the genome, and graphing metagene data.

Dependencies: matplotlib, numpy, pandas, pysam, cPickle
These are all installed in a virtualenv on honey-badger. You can activate it by typing `source /home/susan/susan/bin/activate`. If you are not working in an environment with these dependencies installed, the program will throw an error.

There are 2 ways to run the coverage suite:

**User Script (hand holding version)**  
Type `python user_script_coverage.py` from this folder and it will talk you through the steps you need to take.

**Manual Version**
Type `python coverage_user.py function arguments` from this folder.
Here are the functions available to you through that version.

* **index_and_pileup ( genome_length, bam_files )**
 `python coverage_user.py index_and_pileup 4309328 [folder1/condition1_sorted.bam, folder2/condition2_sorted.bam, folder3/condition3_sorted.bam]`
This generates an index file *[Bamfilename].bai* and a coverage file *coverage.pkl* for each bam file listed.
*You must run this function once per bam file to generate the necessary data files for the other functions*

* **clean_data ( data_folders, control_folder )**
`python coverage_user.py clean_data [folder1, folder2] folder3`
This function subtracts the controlled data from the listed control folder from the data in each of the listed data folders. generating coverage_controlled.pkl files in each of the experiment folders listed, eg folder2/folder2_coverage_controlled.pkl
*If you don't want to subtract any control data, make sure to set out_of_control=True where it is an option. *

* **graph_range ( genome_range, data_folders, out_of_control=False )**
 `python coverage_user.py graph_data_by_range [0,1000] [folder1,folder2]`
This function generates line graphs for the given range of data for the data in each of the files listed.  
eg. coverage_start_end/coverage_controlled_enclosing_start_end.png The new enclosing folder is created if absent.
If the data does not have a control, add `True` at the end:
`python coverage_user.py graph_data_by_range [0,1000] [folder1, folder2] True`
The word "controlled" will be omitted from the filenames of the graphs generated.
*If you are getting a syntax error, make sure there are no spaces between the brackets, ie [0,1000], not [ 0,1000 ]*


* **graph_all(data_folders, xrange=100000, basey=2, ymin=2e5, out_of_control=False)**
`python coverage_user.py graph_all [folder1, folder2]`
Graphs full genome coverage for data in folders listed. You can change the resolution on the graphs by editing xrange. For the whole genome on one graph, xrange should equal the genome length. For many, smaller graphs, xrange should be smaller. Basey is the log scale of the y axis. out_of_control is True if your data does not have a control. Bar graphs are saved as
You may choose to set some but not all optional parameters, eg:
`python coverage_user.py graph_all [folder1,folder2] xrange=40932010 out_of_control=True` 

* **graph_metagene(CDS_file, data_folders, out_of_control=False, _range=1000)**
`python coverage_user.py graph_metagene CDS_list.csv [folder1,folder2]`
Graphs metagene analysis for data in the given data folders, according to the CDS file given (columns must match the sample provided).
A single graph with name folder_metagene_+-range.png will be created for each data folder given.
xrange indicates the number of bases upstream and downstream of the +1 site that are to be included in the graph. This parameter may be changed as follows:
`python coverage_user.py graph_metagene CDS_list.csv [folder1, folder2] _range=2000`

