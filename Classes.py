# -*- coding: utf-8 -*-
# All the Classes

import numpy as np


class dataClass():
    def __init__(self, name, data):
        self.name = name
        self.data = np.array(data)
        self.error = 0.0
        self.fitCoeff = 1.0
        self.color = 'b'
        self.style = 'scatter'
