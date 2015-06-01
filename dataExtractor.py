# -*- coding: utf-8 -*-
# Data extractor for BB

import numpy as np


def extractData(filename, columnNumber=3):
    # Extracts data from a file with a known number of columns of data
    data = np.genfromtxt(filename)
    outputList = []

    for i in data:
        outputList.append(i[columnNumber])  # ouputs only 1 column of data

    return outputList
