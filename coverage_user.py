import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['agg.path.chunksize'] = 10000
import matplotlib.pyplot as plt
import cPickle as pkl
import pandas as pd
import numpy as np
import sys
import pysam
import os
	
def main():
        function=sys.argv[1]
        args=sys.argv[2:]
        for i in range(len(args)):
                if args[i][0]=="[":
                        args[i]=[a.strip() for a in args[i].replace("]","").replace("[","").split(",")]
                        for arg in args[i]:
                                if arg.isdigit():
                                        arg=int(arg)
                        print args[i]
        if len(args)==1:
                eval(function)(args[0])
        elif len(args)==2:
                eval(function)(args[0],args[1])
        elif len(args)==3:
                eval(function)(args[0],args[1],args[2])

def index_and_pileup(genome_len, bam_files):
	for bam_file in bam_files:
		print "processing", bam_file
	        pysam.index(bam_file)
		coverage=[0]*genome_len
	        f=pysam.Samfile(bam_file, 'rb')
        	for pileup_col in f.pileup():
                	coverage[pileup_col.pos]=pileup_col.n
		enclosing_folder="/".join(bam_file.split("/")[:-1])
		immediate=bam_file.split("/")[-2]
        	pkl.dump(coverage,open(enclosing_folder+"/"+immediate+"_coverage.pkl",'wb'))

def clean_data(enclosingfolders, controlfolder):
	immediate=controlfolder.split("/")[-1]
	control=pkl.load(open(controlfolder+"/"+immediate+"_coverage.pkl", 'rb'))
	for enclosing in enclosingfolders:
		immediate=enclosing.split("/")[-1]
		data=pkl.load(open(enclosing+"/"+immediate+"_coverage.pkl", 'rb'))
		for i in range(len(data)):
			data[i]=data[i]-control[i]
		pkl.dump(data, open(enclosing+"/"+immediate+"_coverage_controlled.pkl", "wb"))

def load_clean_data(filename):
	print "loading", filename
	data=pkl.load(open(filename, "rb"))
	return data

def graph_range(genomerange, enclosings, out_of_control=False):	
	for enclosing in enclosings:
		immediate=enclosing.split("/")[-1]
		if out_of_control:
			data=load_clean_data(enclosing+"/"+immediate+"_coverage.pkl")
		else:
			data=load_clean_data(enclosing+"/"+immediate+"_coverage_controlled.pkl")
		start=genomerange[0]
		end=genomerange[1]
		d=os.path.dirname(enclosing)+'/coverage_'+str(start)+"_"+str(end)
		if not os.path.exists(d):
			os.makedirs(d)
		subdata=data[start:end]
		fig, ax=plt.subplots()
		plt.plot(range(len(subdata)), subdata,linestyle="solid")
		ax.set_ylim(ymin=0)
		ax.set_xticks(range(0,len (subdata), 500))
		ax.set_xticklabels( range(start, end, 500))
		if out_of_control:
			plt.savefig(d+'/coverage_'+immediate+"_"+str(start)+"_"+str(end))
		else:
			plt.savefig(d+'/coverage_controlled_'+immediate+"_"+str(start)+"_"+str(end))
		plt.clf()
		plt.close()
		print "Graphed", enclosing, "from", start, "to", end

def graph_all(enclosings, _range=100000, basey=2, ymin=2e5, out_of_control=False):
	for enclosing in enclosings:
		d=os.path.dirname(enclosing)
		immediate=enclosing.split("/")[-1]
		print enclosing, "about to load data"
		if out_of_control:
			cleaned_data=load_clean_data(enclosing+"/"+immediate+"_coverage.pkl")
		else:
			cleaned_data=load_clean_data(enclosing+"/"+immediate+"_coverage_controlled.pkl")
		ymax=max(cleaned_data)
		for i in range(_range, len(cleaned_data), _range):
			start=max(0, i-_range)
			end=min(i,len(cleaned_data))
			fig, ax=plt.subplots()
			ax.set_xticks(range(0,min(_range, (end-start)), _range/3))
			ax.set_xticklabels(range(start, end, _range/3))
			ax.set_yscale('log', basey=basey)
			ax.set_ylim([ymin,ymax])
			a=plt.plot(cleaned_data[start:end])
			ax.fill_between(a[0].get_xdata(), 0,cleaned_data[start:end])
			if out_of_control:
				plt.savefig(enclosing+"/"+immediate+"_"+str(start)+"_"+str(end)+"_coverage")
			else:
				plt.savefig(enclosing+"/"+immediate+"_"+str(start)+"_"+str(end)+"_controlled_coverage")
			print "graphed " , enclosing," ",str(start),"-",str(end)
			plt.close()

def parse_plus1s(filename):
	df=pd.read_csv(filename)
	plus1s=df["Left-End-Position"].where(df["Transcription-Direction"]=="+", other=df["Right-End-Position"])
	return plus1s

def plot_plus1s(plus1s, data, enclosing, _range=1000):
	aligneddata=[0]*((2*_range)+1)
	counts=[0]*((2*_range)+1)
	for index in plus1s:
		try:	
			index=int(index)
			for i in range(-_range,_range+1):
				aligneddata[i+_range]+=data[(i+index)%len(data)]
				counts[i+_range]+=1
		except ValueError:
			pass 
	averageddata=[]
	for i in range(len(aligneddata)):
		averageddata.append(aligneddata[i]/float(counts[i]))
	fig, ax=plt.subplots()
	ax.plot(range(len(averageddata)),averageddata,linestyle="solid")
	plt.xticks(np.arange(0.0, float(_range), _range/5),np.arange(-_range,_range,_range/5))
	d=os.path.dirname(enclosing)
	immediate=enclosing.split("/")[-1]
	plt.savefig(d+"/"+immediate+'_metagene+-'+str(_range))
	print "plotted metagene analysis for ",enclosing

def graph_metagene(plus1sfile, folders, out_of_control=False, _range=1000):
	for enclosing in folders:
		immediate=enclosing.split("/")[-1]
		if out_of_control:
	                cleaned_data=load_clean_data(enclosing+"/"+immediate+"_coverage.pkl")
		else:
			cleaned_data=load_clean_data(enclosing+"/"+immediate+"_coverage_controlled.pkl")
		print "loaded data", enclosing
		plus1s=parse_plus1s(plus1sfile)
		plot_plus1s(plus1s,cleaned_data, enclosing, _range)

