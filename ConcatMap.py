#!/usr/bin/env python

"""
Copyright 2019  Daryl Gohl

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import interp1d
import numpy as np
import pysam
import os
import argparse
__version__ = "1.0"

def get_args(x):
    x.add_argument("-i", "--input_file",
                   type = str,
                   default = None,
                   metavar = '',
                   help="Input path for sam file [required].")
    x.add_argument("-r", "--reference_length",
                   type = int,
                   default = None,
                   metavar = '',
                   help="Input reference length [required].")
    x.add_argument("-m", "--min_length",
                   type = int,
                   default = 100,
                   metavar = '',
                   help="Minimum mapped read length to plot (default: 100)")
    x.add_argument("-o", "--output_dir",
                   type = str,
                   default = '',
                   metavar = '',
                   help = "Output directory for plot (default: same folder as input file)") 
    x.add_argument("-l", "--line_spacing",
                   type = float,
                   default = 0.02,
                   metavar = '',
                   help = "Radial spacing of each read on plot (default 0.2)")
    x.add_argument("-w", "--line_width",
                   type = float,
                   default = 0.75,
                   metavar = '',
                   help = "Line width of each read on plot (default 0.75)")
    x.add_argument("-c", "--circle_size",
                   type = float,
                   default = 0.45,
                   metavar = '',
                   help = "Size of central circle (default 0.45)")
    x.add_argument("-s", "--fig_size",
                   type = float,
                   default = 10,
                   metavar = '',
                   help = "Size of figure (default 10)")
    x.add_argument("-f", "--figure_format",
                   type = str,
                   default = "pdf",
                   metavar = '',
                   help = "Format of saved figure, supported formats: eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff.")

#####################    
    args = x.parse_args()
    return args

# Parses command line arguments
argparser = argparse.ArgumentParser(description = "Map and plot reads against a circular reference (v" + __version__ +")\n" + \
                                              "by Daryl Gohl\n" + \
                                              "This program takes in SAM files of sequencing reads mapped to a concatenated reference sequence (two copies of the reference sequence repeated in tandem) and the reference sequence length (length of the original presumed circular reference sequence) and outputs a plot of the mapped reads against a circularized reference.",
                                add_help = True, 
                                epilog ='')
args = get_args(argparser)
args = vars(args)

# possible input args
filename = args['input_file']
fname = os.path.split(filename)[1]
ref_length = args['reference_length']
incrementor = args['line_spacing']
width = args['line_width']
out_folder = args['output_dir']
Min_length = int(args['min_length'])
Fig_size = args['fig_size']
start_coord = args['circle_size']
Fig_format = args['figure_format']

if out_folder == '':
  out_dir = os.path.dirname(filename)
else:
  out_dir = out_folder

#Set output file name
file = fname[:-4] + "." + Fig_format
fig_name = os.path.join(out_dir,file)

#SAM File
samfile = pysam.AlignmentFile(filename, "r")

#reads
reference_start = []
reference_end = []
for r in samfile.fetch(until_eof=True):
    x = r.query_alignment_start #position of start of alignment in read
    y = r.query_alignment_end #position of end of alignment in read
    z = r.query_alignment_length #length of alignment in read
    x1 = r.reference_start #position of start of alignment in reference
    y1 = r.reference_end #position of end of alignment in reference
    z1 = r.reference_length #length of start of alignment in reference
    if z1 > Min_length:
        reference_start.append(x1)
        reference_end.append(y1)

#Dedup reference
break_span = []
ref_start_collapsed = [] 
ref_end_collapsed = []
for i, item in enumerate(reference_start):
    if item < ref_length and reference_end[i] > ref_length:
        break_span.append(True)
    else:
        break_span.append(False)
    if item > ref_length:
        ref_start_collapsed.append(item-ref_length)
    else:
        ref_start_collapsed.append(item)
    if reference_end[i] > ref_length:
        ref_end_collapsed.append(reference_end[i]-ref_length)
    else:
        ref_end_collapsed.append(reference_end[i])

#Convert to polar coordinated (0-360)
deg_min = 0
deg_max = 360
deg_per_base = 360.0/ref_length
deg_start = []
deg_end = []
for i, item in enumerate(ref_start_collapsed):
    deg_start.append(item*deg_per_base)
    deg_end.append(ref_end_collapsed[i]*deg_per_base)


#Plotting data
with plt.style.context("seaborn-white"):
    fig = plt.figure(figsize=(Fig_size,Fig_size))
    #my_dpi = 300
    #fig = plt.figure(figsize=(2400/my_dpi, 2400/my_dpi), dpi=my_dpi)
    ax = fig.add_subplot(111, projection="polar")
    ax.grid(False)
    ax.set_rticks([])
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_theta_zero_location('N')
    ax.set_facecolor('white')
    ax.axis('off')

    # Connect two points with a curve
    for curve in [[[0, 360], [(start_coord-2*incrementor), (start_coord-2*incrementor)]]]:
        curve[0] = np.deg2rad(curve[0])
        x = np.linspace( curve[0][0], curve[0][1], 500)
        y = interp1d( curve[0], curve[1])( x)
        ax.plot(x, y, linewidth = 5)
    for curve in [[[0, 0], [0.0, 0.0]]]:
        curve[0] = np.deg2rad(curve[0])
        x = np.linspace( curve[0][0], curve[0][1], 500)
        y = interp1d( curve[0], curve[1])( x)
        ax.plot(x, y)

    #reads    
    for i, item in enumerate(break_span):
        if item == False:
            for curve in [[[deg_start[i], deg_end[i]], [start_coord, start_coord]]]:
                curve[0] = np.deg2rad(curve[0])
                x = np.linspace( curve[0][0], curve[0][1], 500)
                y = interp1d( curve[0], curve[1])( x)
                ax.plot(x, y, color='grey', linewidth = width)
        else:
            for curve in [[[360, deg_start[i]], [start_coord, start_coord]]]:
                curve[0] = np.deg2rad(curve[0])
                x = np.linspace( curve[0][0], curve[0][1], 500)
                y = interp1d( curve[0], curve[1])( x)
                ax.plot(x, y, color='grey', linewidth = width)
            for curve in [[[deg_end[i], 0], [start_coord, start_coord]]]:
                curve[0] = np.deg2rad(curve[0])
                x = np.linspace( curve[0][0], curve[0][1], 500)
                y = interp1d( curve[0], curve[1])( x)
                ax.plot(x, y, color='grey', linewidth = width)
        start_coord = start_coord + incrementor
    
plt.savefig(fig_name,bbox_inches='tight')
