import numpy
import pandas as pd
import scipy
import os
import sys
from ibidas import *
import yaml

result = {}
autofit = None

files = sys.argv[1:4]
result_file = sys.argv[5]
result_file2 = sys.argv[6]

for fname, colname in zip(files, ['auto', 'chry','chrx','chrm']):
    #read data
    res = pd.read_csv(fname, delimiter='\t')
    kmers = res.iloc[:,0]

    #gc counts
    gc = []
    for kmer in kmers:
        g = len([e for e in kmer if e in ['C','G']])
        gc.append(g - 16)
    gc = numpy.array(gc)


    counts = res.iloc[:,1].to_numpy()
    
    #correction for GC
    if colname == 'auto':
        autofit = scipy.stats.linregress(gc, numpy.log(counts + 0.01))
        result['slope'] = float(autofit.slope)
        result['intercept'] = float(autofit.intercept)
        result['rvalue'] = float(autofit.rvalue)
        print(autofit)
    elif colname == 'chry':
        counts = counts - 1  #correction for earlier processing step
    correction = autofit.slope * gc

    #calculate median and percentiles
    result[colname] = float(numpy.median(counts))
    #GC correction does not seem to add much, disabled for now.
    #counts = numpy.exp(numpy.log(counts + 0.01) - correction) - 0.01
    #counts[counts < 0] = 0.001
    result['dist_' + colname] = [float(e) for e in numpy.percentile(counts, numpy.linspace(0, 100, 101))]


result['xratio'] = float(2.0 * result['chrx'] / result['auto'])
result['yratio'] = float(2.0 * result['chry'] / result['auto'])

result['yyratio'] = float(2.0 * result['chry'] / (result['chrx'] + result['chry']))

result['dist_chrx_ratio'] = [float(e) for e in (numpy.array(result['dist_chrx']) / numpy.array(result['dist_auto']))]
result['dist_chry_ratio'] = [float(e) for e in (numpy.array(result['dist_chry']) / numpy.array(result['dist_auto']))]



result['sex'] = 'F' if result['yyratio'] < 0.05 else 'M'
Save(result,result_file)
with open(result_file2, 'w') as f:
    yaml.dump(result, f, default_flow_style=False)
