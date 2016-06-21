from coverage_user import *
import os

def main():
	yeses=['y','Yes','Y','yes']
	nos=['n', 'no', 'N', 'No']
	filenames=[]
	firsttime=raw_input("Hello! Is this your first time running this program for this dataset? y/n\n")
	while firsttime not in yeses+nos:	
		firsttime=raw_input("please type a valid response, y or n")
	if firsttime in yeses:
		filenames=raw_input("Yes? Let's create some index files and pileups of coverage. Type a comma separated list of your sorted bam files (including paths), then press enter.\n")
		genome_len=int(raw_input("How long is the genome?\n"))
		print "This will take about 5 minutes for each file"
		filenameslist=[name.strip(" ") for name in filenames.split(",")]

		index_and_pileup(genome_len,filenameslist)

	controls=raw_input("Ok. Have you already subtracted any controls you want to? Type ! if you will not be subtracting any controls. Type y/n/!\n")
	while controls not in yeses+nos+['!']:
		controls=raw_input("please type a valid response, y, n, !\n")
	if controls in nos:
		controlfile=raw_input("What folder contains your control file?\n")
		filenames=raw_input("What folders contain the data you want to subtract this control from? Please type a comma separated list.\n")	
		filenameslist=[name.strip(" ") for name in filenames.split(",")]
		print "this will take just a minute"
		clean_data(filenameslist, controlfile)
	if controls=="!":
		out_of_control=True
	else:
		out_of_control=False

	a=raw_input("Do you want to graph coverage for a particular region of the genome? y/n\n")
	while a not in yeses+nos:
		a=raw_input("Please type a valid response, y or n\n")
	if a in yeses:
		if filenames:	
			print "You listed these filenames: ", filenameslist
			b=raw_input("Do you want to use these again? y/n\n")
			while b not in yeses+nos:
				b=raw_input("please type a valid response, y, n\n")
			if b in nos:
				filenames=raw_input("What folders contain your data?\n")
			filenameslist=[name.strip(" ") for name in filenames.split(",")]
		else:
			filenames=raw_input("What folders contain your data?\n")
                        filenameslist=[name.strip(" ") for name in filenames.split(",")]
		still_going=True
		while still_going:
			start=int(raw_input("What index do you want to start at?\n"))
			end=int(raw_input("What index do you want to end at?\n"))
			graph_range([start,end], filenameslist, out_of_control)      	
			c=raw_input("Do you want to graph other ranges for the same files? y/n\n")
			while c not in yeses+nos:
				c=raw_input("please type a valid response, y, n\n")
			if c in yeses:
				still_going=True	
			else:
				still_going=False
	a=raw_input("Do you want to graph whole genome coverage? y/n\n")
	while a not in yeses+nos:
                a=raw_input("Please type a valid response, y or n\n")
        if a in yeses:	
		if filenames:
                        print "You listed these filenames: ", filenameslist
                        b=raw_input("Do you want to use these again? y/n\n")                        
			while b not in yeses+nos:
                                b=raw_input("please type a valid response, y, n\n")
                        if b in nos:
                                filenames=raw_input("What folders contain your data?\n")
                         	filenameslist=[name.strip(" ") for name in filenames.split(",")]
                else:   
                        filenames=raw_input("What folders contain your data?\n")
                        filenameslist=[name.strip(" ") for name in filenames.split(",")]
		graph_all(filenameslist, out_of_control=out_of_control)
	
	a=raw_input("Do you want to graph metagene analysis? y/n \n")
	while a not in yeses+nos:
                a=raw_input("Please type a valid response, y or n\n")
        if a in yeses:
		CDS_list=raw_input("What is the filename of your cds list? It must be a csv and formatted like the example.\n")
		if filenames:
                        print "You listed these filenames: ", filenameslist
                        b=raw_input("Do you want to use these again? y/n\n")
			while b not in yeses+nos:
                                b=raw_input("please type a valid response, y, n\n")
                        if b in nos:
                                filenames=raw_input("What folders contain your data?\n")
	                        filenameslist=[name.strip(" ") for name in filenames.split(",")]
		else:
                        filenames=raw_input("What folders contain your data?\n")
		        filenameslist=[name.strip(" ") for name in filenames.split(",")]
		graph_metagene(CDS_list, filenameslist, out_of_control=out_of_control)
	print "You're all done. Goodbye"


main()
