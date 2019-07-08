# ConcatMap
This program takes in SAM files of sequencing reads mapped to a concatenated reference sequence (two copies of the reference sequence repeated in tandem) and the reference sequence length (length of the original presumed circular reference sequence) and outputs a plot of the mapped reads against a circularized reference.

## Prerequisites
Python 2.7
matplotlib
scipy
numpy
pysam

## Usage
ConcatMap.py [-h] [-i] [-r] [-m] [-o] [-l] [-w] [-c] [-s] [-f]

optional arguments:

    -h, --help            show this help message and exit
  
    -i , --input_file     Input path for sam file [required].
  
    -r , --reference_length 
                        Input reference length [required].
  
    -m , --min_length     Minimum mapped read length to plot (default: 100)
  
    -o , --output_dir     Output directory for plot (default: same folder as
                        input file)
  
    -l , --line_spacing   Radial spacing of each read on plot (default 0.2)
  
    -w , --line_width     Line width of each read on plot (default 0.75)
  
    -c , --circle_size    Size of central circle (default 0.45)
  
    -s , --fig_size       Size of figure (default 10)
  
    -f , --figure_format 
                        Format of saved figure, supported formats: eps, jpeg,
                        jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif,
                        tiff.

## Usage example
ConcatMap.py -i <PathToFile/InputFileName> -r <ReferenceLength> 
