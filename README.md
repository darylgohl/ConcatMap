# ConcatMap
This program takes in a fasta reference file and a fastq sequencing file and concatenates the reference sequence (two copies of the reference sequence repeated in tandem) and then uses minimap2 to map reads to the concatenated reference, generating a .sam file that is used for subsequent plotting.

## Prerequisites
This program relies on minimap2, an external tool for sequence alignment. You must install it separately before using ConcatMap. Follow the instructions on the minimap2 [https://github.com/lh3/minimap2] GitHub page to install it for your operating system.

## Install
You can install ConcatMap directly from the GitHub repository using pip. It's recommended to do this within a virtual environment.
pip install git+https://github.com/darylgohl/ConcatMap.git

## Usage
ConcatMap [-h] [-q] [-r] [-o] [-n] [-m] [-u] [-l] [-w] [-c] [-s] [-x] [-f]

Map and plot reads against a circular reference (v1.2) by Daryl Gohl This
program takes in SAM files of sequencing reads mapped to a concatenated
reference sequence (two copies of the reference sequence repeated in tandem)
and the reference sequence length (length of the original presumed circular
reference sequence) and outputs a plot of the mapped reads against a
circularized reference.

optional arguments:

  -h, --help            show this help message and exit

  -q , --fastq_file     Input path for fastq sequencing reads file [required].

  -r , --fasta_file     Input path for fasta reference file [required].

  -o , --output_dir     Output directory for sam file and plot (default: same folder as input fastq file)

  -n , --output_file    Output name for sam file and plot (default: same name as input fastq file)

  -m , --min_length     Minimum mapped read length to plot (default: 100)
  
  -u , --unsorted     Plot from unsorted sam file (default: False - plot sorted sam file)

  -l , --line_spacing   Radial spacing of each read on plot (default 0.2)

  -w , --line_width     Line width of each read on plot (default 0.75)

  -c , --circle_size    Size of central circle (default 0.45)

  -s , --fig_size       Size of figure (default 10)

  -x , --clip           Plot clipped portion of reads (default False)

  -f , --figure_format	Format of saved figure, supported formats: eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff.
                        
## Usage example
ConcatMap -q <PathToFASTQFile/InputFileName> -r <PathToFASTAReferenceFile/ReferenceFileName> -m 10000 -f "png"

