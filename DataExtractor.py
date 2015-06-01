# -*- coding: utf-8 -*-
# Data extractor for BB
import numpy as np

def extractData(filename, columnNumber = 3):
    data = np.genfromtxt(filename)
    outputList = []
    
    for i in data:
        outputList.append(i[3])
        
    return outputList
