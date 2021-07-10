DESCRIPTION:
   DiSCo: Dsr-dependent dIssimilatory Sulfur metabolism Classification tOol
	- Automatic detection and prediction of enzyme type of proteins involved in Dsr-dependent dissimilatory sulfur metabolism.


All programs are written in perl (version 5) and dependent on HMMER3. 
The script requires hmmsearch to be in your PATH. 

This tool runs a set of FASTA formatted protein sequences against the curated set of hidden Markov models (HMMs).
Please report bugs and requests to: towhom?@domain.com


1) Download DiSCo:
	https://github.com/Genome-Evolution-and-Ecology-Group-GEEG/DiSCo

2) Get HMMER3
	Download HMMER3 from http://hmmer.org/ and follow its instructions.

3) Run DiSCo to detect proteins and predict the enzyme type:
	perl DiSCo.pl -i example.faa 

4) Output format, column description:
Sequence identifier
HMM name
Predicted enzyme type
Input sequence length
HMM length
Score
Bias
E-value
Accuracy
HMM alignment length 
Conditional E-value
Independent E-value
Domain score
Domain bias
Current domain
Total number of domains
Domain start with respect to the input sequence
Domain end with respect to the input sequence
Input sequence annotation
Notes: DsrC-type (1 or 2 cysteines)

############################################################################
USAGE:
  ./DiSCo.pl -i <filename> [-p <file_prefix>] [-o]
  [-d </path/to/dir/>]
  [-m </path/filename>] [-f </path/filename>]
  [-t </path/filename>]
  [-s <number>] [-n <number>] 
  [-q] [-h] [-v]

 example:
  perl DiSCo.pl -i example.faa  (runs DiSCo with default settings)
  perl DiSCo.pl -i example.faa -p example  (runs DiSCo with given filename prefix)
  perl DiSCo.pl -i example.faa -p example -s 1 -o (runs DiSCo with given filename prefix, specifies outfile separator, and produce sequence outfile)

Standard settings:
  -i <filename>		input protein sequence file in FASTA format

Additional settings:
  -p <prefix>		filtered DiSCo hitlist - filename prefix
  -s <number>		separator for hitlist
   0 = TAB=\t		[default: tab separated]
   1 = ,
   2 = ;
   3 = SPACE=" "
	
  -d </path/to/dir>			absolute path to outfile directory	[optional, default current dir: ./]
  -o 						DiSCo hits - sequences output in FASTA format			[optional]
  -n <number>				number of threads for hmmsearch		[optional, default = 1]
  -m 						alternative location and filename of DiSCo HMM library	[optional, default current dir: ./]
  -f 						alternative location and filename of filter script	[optional, default current dir: ./]
  -t </path/thresholds.txt>	tab seperated table with user-defined cut-offs: column 1: model, column 2: score, column 3: E-value

General settings:
  -h 	show help and usage
  -q 	quiet, non verbose
  -v 	version


