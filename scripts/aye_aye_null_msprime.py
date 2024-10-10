import sys
import pandas as pd
import math
import os
import numpy as np
import msprime
import tskit
import libsequence
import allel
import argparse


#Example input
#python aye_aye_null_msprime.py -region 1 -seq_len 50000 -num_replicates 100 \
#-outPath "/home/vivak/"

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about length of region and sample size')
parser.add_argument('-region', dest = 'region', action='store', nargs = 1, type = int, help = 'Region')
parser.add_argument('-num_replicates', dest = 'num_replicates', action='store', nargs = 1, type = int, help = 'Number of replicates')
parser.add_argument('-outPath', dest = 'outPath', action='store', nargs = 1, type = str, help = 'path to output files (suffixes will be added)')

args = parser.parse_args()
region = args.region[0]
num_replicates = args.num_replicates[0]
outPath = args.outPath[0]


chr_lens = { 1:18602764, 2:16060938, 3:17214299, 4:14784439, 5:12151006, 6:12101961, 7:8717938, 8:9958554, 10:6331343, 11:2950556, 12:2046229, 13:1970739, 14:1919907, 15:1563188 }

N = 23389/2
t_bot = 1133/2
N_north_bot =  2691/2 
N_other_bot =  6585/2

t_decline = 7
N_north_current = 1050/2
N_other_current = 2570/2

r_north = (N_north_current/N_north_bot)**(1/t_decline)-1.008641
r_other = (N_other_current/N_other_bot)**(1/t_decline)-1.008641

def aye_aye_demog(seq_len, num_replicates):
    demography = msprime.Demography()
    demography.add_population(
        name="aye_aye",
        description="Other Aye-Aye population",
        initial_size=N_other_current,
        growth_rate=r_other,
        default_sampling_time=0, 
        initially_active=True,
    )
    demography.add_population(
        name="aye_aye_anc",
        description="Ancestral Aye-Aye population",
        initial_size=N,
        growth_rate=0,
    )
    #Add events
    demography.add_population_parameters_change(time=t_decline, growth_rate=0, population="aye_aye")
    demography.add_population_split(time=t_bot, derived=["aye_aye"], ancestral="aye_aye_anc")
    demography.sort_events()

    ancestry_reps = msprime.sim_ancestry(
        {"aye_aye": 5}, 
        demography=demography, 
        sequence_length = seq_len,
        recombination_rate = 1e-8,
        num_replicates=num_replicates)

    for ts in ancestry_reps:
        mutated_ts = msprime.sim_mutations(ts, rate=1.52e-8)
        yield mutated_ts


ts_list = []
for replicate_index, ts in enumerate(aye_aye_demog(chr_lens[region], num_replicates)):
    ts_list.append(ts)


for rep, ts in enumerate(ts_list):
    lst = []
    #Loop through variants
    for var in ts.variants():
        #Obtain counts of ancestral and derived alleles
        unique, counts = np.unique(var.genotypes, return_counts=True)
        #Append list of positions and MACs to main list
        lst.append([int(var.site.position), counts.min()])
    #Convert to df     
    df = pd.DataFrame(lst, columns=['position', 'x'])
    df['n'] = 10
    df['folded'] = 1
    df = df[df.x<10]

    df2 = pd.DataFrame(lst, columns=['physPos', 'x'])
    df2['n'] = 10
    df2['genPos'] = 'NA'
    df2 = df2[['physPos', 'genPos', 'x', 'n']]
    df2 = df2[df2.x<10]
    
    df.to_csv(outPath + "/scaffold" + str(region) + "_rep" + str(rep) + "_SF2.aff", header=True, index=False, sep='\t')
    df[['position']].to_csv(outPath + "/scaffold" + str(region) + "_rep" + str(rep) + "_SF2.grid", header=False, index=False, sep='\t')
    df2.to_csv(outPath + "/scaffold" + str(region) + "_rep" + str(rep) + "_BM.aff", header=True, index=False, sep='\t')
