#!/usr/bin/env python

# python module to extract longitudes, latitudes and distances, closest to n deltad.
# Gareth Roberts, Imperial College London, 5 Dec 2019

import numpy as np
import pandas as pd

deltad = 2

df = pd.read_csv('distances', delim_whitespace=True, names=['lon', 'lat', 'distance'])
dlen = (np.int(df['distance'].max()))

with open('distances.out', 'a') as f:
    for i in range(0,dlen,deltad): 
        value = np.real(i)
        
        newdf = df.iloc[(df['distance']-value).abs().argsort()[0:1]]    
        newdf.to_csv(f, sep='\t', header=False)
