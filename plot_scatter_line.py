#!/usr/bin/python

import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt
import prettyplotlib as ppl
import seaborn as sns
#from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import *
from pylab import *
#from BrewerColors import *




def main(argv):

    # get options passed at command line
    try:
        opts, args = getopt.getopt(argv, "p:o:t:f:l:x:")
    except getopt.GetoptError:
        #print helpString
        sys.exit(2)


    #print opts

    for opt, arg in opts:
        if opt == '-p':
            params = arg
        elif opt == '-t':
            title = arg
        elif opt == '-o':
            outfile_name = arg
        elif opt == '-f':
            file_names = arg
        elif opt == '-l':
            labels = arg
        elif opt == '-x':
            x_vals = arg



    file_names = file_names.split(",")
    #print file_names
    labels = labels.split(",")
    #print x_vals
    x_vals = x_vals.split(",")
    #print x_vals
    x_vals = [int(x) for x in x_vals]

    #colors = sns.color_palette("husl", 20)

    colour_dict = {'DiamondMegan_filtered':(190,190,190),
             'CLARK':(178,223,138),
             'CLARK-S':(51,160,44),
             'DiamondEnsemble':(227,26,28),
             'BlastEnsemble':(251,154,153),
             'Kraken':(255,127,0), 
             'Kraken_filtered':(253,191,111),
             'PhyloSift':(102,205,0),
             'PhyloSift_filtered':(127,255,0),
             'MetaPhlAn':(148,0,211),
             'LMAT':(176,48,96),
             'GOTTCHA':(106,61,154),
             'MetaFlow':(177,89,40),
             'Community':(0,0,0),
             'DiamondMegan_filtered+Kraken_filtered':(202,178,214),
             'CLARK+GOTTCHA':(255,233,0),
             'BlastMegan_filtered+LMAT':(33,160,160),
             'BlastMegan_filtered':(166,206,227),
             'BlastMegan_filtered_liberal':(31,120,180),
             'NBC':(255,0,255)}
    for tool in colour_dict:
        colour_dict[tool] = tuple(map(lambda x: x/255., colour_dict[tool]))

    #colors_255=[(102,205,0),(176,48,96),(51,160,44),(255,127,0),(178,223,138),(253,191,111),(255,0,255),(127,255,0),(227,26,28),(31,120,180),(251,154,153),(106,61,154),(148,0,211),(202,178,214),(177,89,40),(166,206,227),(190,190,190)]

    #colors_1 = []
    #for (r,g,b) in colors_255:
    #    colors_1.append((float(r)/255.0, float(g)/255.0, float(b)/255.0))

    #print colors_1

    #colors = colors_1

    figure = plt.figure(figsize=(9, 5), dpi=300)    
    font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 8}

    rc('font', **font)  # pass in the font dict as kwargs


    ##### line plot ######

    lineplt = figure.add_subplot(121)

    data = []
    colors = []
    for file_name,label in zip(file_names,labels):
        colors.append(colour_dict[label])
        data_points = list(np.genfromtxt(file_name, usecols=0))
        print label,data_points
        data.append(data_points)
    data_max = []
    for data_points in data:
        data_max.append(int(max(data_points)))

    data_percent = []
    for i, data_points in enumerate(data):
        percents = []
        for val in data_points:
            per = float(val) / data_max[i] * 100
            percents.append(per)
        data_percent.append(percents)

    print data_max

    print labels

    for i, data_percents in enumerate(data_percent):
        ppl.plot(x_vals[0:len(data_percents)], data_percents, '-', color=colors[i])

        

    xticks(x_vals, rotation=60)
    lineplt.set_xlabel("Input dataset size (Million Reads)")
    yticks(np.arange(0,105,10))
    ylim(top=110)
    lineplt.set_ylabel("Percent max species recovered")
    #plt.show()


    barplt = figure.add_subplot(122)


    labels, data_max, colors = [list(x) for x in zip(*sorted(zip(labels, data_max, colors), key=lambda pair: pair[1]))]

    ppl.barh(np.arange(0,len(labels)), data_max, log=1, color=colors)
    yticks(np.arange(0,len(labels)), labels)
    barplt.yaxis.label.set_verticalalignment('top')
    barplt.set_xlabel("Maximum number of species predicted")

    #plt.legend(labels, bbox_to_anchor=(0.98,0.68))
    plt.tight_layout()

    plt.savefig(outfile_name, format='pdf')

if __name__ == "__main__":
    main(sys.argv[1:])



















