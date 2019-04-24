'''
This script takes the results of RepeatExplorer (assembly/contigs and assembly/contigs.prof) and gets a list of contigs from the RepeatExplorer analysis that are highly repetitive.

Use to screen the file contigs.prof for sequences with large numbers of repeats.

It outputs a file called Output Files:
	Best_Contigs.txt - a list of contigs and a number representing the depth of reads in that contig, as a proxy for the repetitiveness of that contig.  The cutoff for the depth needed is 1000, which can be changed inside the script.
	PrimerTargets.txt- a file that contains, for each of the contigs returned in Best_Contigs.txt, the original contig sequence, the depth profile for the entire sequence, the depth profile for the highly repetitive region and the sequence of the highly repetitive region.  The depth considered to be highly repetitive is by default 1000 or the maximum depth +/- 25% of the maximum depth, whichever is less.  This can be changed within the script.  See RepeatExplorer documentation for details of the profile files.  

Usage:
	Run from wherever the contigs and contigs.prof file is (e.g. /seqClust/assembly/)
	<python FindPrimerTargets.py>

The sequence of the highly repetitive region from PrimerTargets.txt can be used as the template for primer design in Primer3 or any good primer design software.

There is a function in the script called \blast() that will blast the repetative regions. It comes commented out becasue it is a little time consuming, but you can uncomment it. For this to work, you must have Bioython installed.

RepeatExplorer information is found here: http://repeatexplorer.org/

The principles of RepeatExplorer approach are described in:
Novak, P., Neumann, P., Macas, J. (2010) - Graph-based clustering and characterization of repetitive sequences in next-generation sequencing data. BMC Bioinformatics 11: 378.
Novak, P., Neumann, P., Pech, J., Steinhaisl, J., Macas, J. (2013) - RepeatExplorer: a Galaxy-based web server for genome-wide characterization of eukaryotic repetitive elements from next generation sequence reads. Bioinformatics

Repeat Explorer uses RepeatMasker and Repbase:
Jurka, J., Kapitonov, V.V., Pavlicek, A., Klonowski, P., Kohany, O., Walichiewicz, J. (2005) Repbase Update, a database of eukaryotic repetitive elements. Cytogentic and Genome Research 110:462-467

Conserved Domain Database:
Geer L.Y., Geer R.C., Gonzales N.R., Gwadz M., Hurwitz D.I., Jackson J.D., Ke Z., Lanczycki C.J., Lu F., Marchler G.H., Mullokandov M., Omelchenko M.V., Robertson C.L., Song J.S., Thanki N., Yamashita R.A., Zhang D., Zhang N., Zheng C., Bryant S.H.(2011) CDD: a Conserved Domain Database for the functional annotation of proteins. Nucleic Acids Res. 39(Database issue):D225-9.

Clustering is performed using Louvain method:
Blondel V.D., Guillaume J., Lambiotte R., Lefebvre E. (2008) Fast unfolding of communities in large networks. J. Stat. Mech.: P10008

author: jgrant@smith.edu
date: 20150312

'''
import os
cutoff = 1000 #can change for different low range for repeat counts
def get_high_counts():
	infile = open('contigs.prof','r')
	di = {}
	li = []
	testdi = {}
	numdi = {}
	for line in infile:
		if line[0] == '>':
			seq = line.split()[0].strip('>')
			numdi[seq] = []
		else:
			for num in line.split():
				numdi[seq].append(int(num))
			if  max(numdi[seq]) > cutoff:
				numdi[seq] = max(numdi[seq])
			else:
				 del numdi[seq]

			
	for seq in numdi:
		li.append(int(numdi[seq]))
		try:
			di[int(numdi[seq])].append(seq)
		except:
			di[int(numdi[seq])] = [seq]

	li.sort(reverse = True)
	outfile = open('Best_Contigs.txt','w')
	outfile.write('Maximum Depth,Sequence Name\n')
	outfile.close()
	for num in li:
		if num > cutoff:
			outfile = open('Best_Contigs.txt','a')
			outfile.write(str(num) + ',' + str((di[num][0])) + '\n')
			outfile.close()
	#print (di,li)
	return di,li
			
def make_data_file(di,li):
	infile1 = open('contigs','r')
	infile2 = open('contigs.prof','r')
	seqdi = {}
	seqname = ''
	seq = ''
	for line in infile1:
		if line[0] == '>':
			seqdi[seqname] = seq
			seq = ''
			seqname = line.strip('>').strip()
		else:
			seq+= line.strip()


	profdi = {}
	seqname = ''
	prof = ''
	for line in infile2:
		if line[0] == '>':
			seqname = line.strip('>').split()[0]
		else:
			profdi[seqname] = line.strip()
	
	get_region(di, li, seqdi, profdi)
		


def get_region(di,li, seqdi, profdi):
	outfile = open('PrimerTargets.txt','w')
	for double in li:
		outfile = open('PrimerTargets.txt','a')
		co = .25*int(double) # another place you could adjust to refine results
		seqname = di[double][0]
		seq = seqdi[seqname]
		prof = profdi[seqname].split()
		indx = prof.index(str(double))
		mini = min(1000, int(prof[indx]) - co)
		indexmax = indx
		indexmin = indx
		#print(indexmin, indexmax)
		for count in range(indx, len(prof)):
			if int(prof[count]) >= mini:
				indexmax += 1
			else:
				break
		for count in range(indx, 0 , -1):
			if int(prof[count]) >= mini:
				indexmin -= 1
			else:
				break
		#print(indexmin, indexmax)
		region = seq[indexmin:indexmax]
		prof_region = ' '.join(prof[indexmin:indexmax])
		#print('cutoff = ', mini)
		#print(indx, indexmin, indexmax)
		#print(indx, indexmin, indexmax)
		
		#print (prof[indexmin], prof[indexmax])
		#print(prof_region)
		#print(region)
		outfile.write(seqname + ', max depth: ' + str(double) + ',length: ' + str(len(region)) +   '\n')
		outfile.write('full sequence: ' + seq + '\n')
		outfile.write('full depth profile: ' + profdi[seqname] + '\n')
		outfile.write('deep region profile: ' + prof_region + '\n')
		outfile.write('deep region sequence: ' +region + '\n'+ '\n')
		outfile.close()
		#blast (region,seqname)
	
	
def blast(seq,name):
	try:
		from Bio.Blast import NCBIWWW
		result_handle = NCBIWWW.qblast("blastn", "nt", seq)
		blast_file = open(name + "_seq_blast_result.xml", "w")
		blast_file.write(result_handle.read())
		blast_file.close()
		result_handle.close()	
		print ('Your blast result is in the file "' + name + '_seq_blast_result.xml"  ')
	except:
		print ("if you had biopython installed, I would automatically blast your sequence of interest. ")

		
		
if __name__ == '__main__':
	di,li = get_high_counts()
	make_data_file(di, li)