# RepeatTargetScripts
Scripts to help parse RepeatExplorer output files for use in developing repeat-based diagnostic assays
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
date: 201904212

