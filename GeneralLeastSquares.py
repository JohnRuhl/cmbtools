import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

#------------------------------------------
# Creating Test Data

# making a test xList and errorYList
xList = [x for x in range(10)]

#errorYList = [np.random.normal(0, 0.1) for y in range(10)]
errorYList = [1.0 for i in range(10)]

# creating test Ylists. These arrays will actually be generated from given functions Y1, Y2, Y3, Y4, etc.
Y1List = [x**1.9 for x in xList]
Y2List = [5.67**(10**-8)*(x**4) for x in xList]
#Y3List = [5,6,7,8,9,10,1,2,3,4]
#Y4List = [4,3,2,1,10,9,8,7,6,6]

yMeasured = [y*np.random.normal(1,0.1) for y in Y1List]


#------------------------------------------
# Matrix Function


# this just makes it easier to refer to all the YLists
YList = [Y1List, Y2List] #[Y1List, Y2List, Y3List, Y4List]

# Creates the A matrix for use in determining the constants
def matrixFunction(functionList, errorYList):
    columns = len(functionList)
    rows = len(errorYList)
    # initialize the matrix with float values of 0
    resultMatrix = np.matrix([[0.0 for i in functionList] for i in errorYList])
    # initialize the list used as temporary storage for the row values
    tempList = np.empty(columns)
    
    for i in range(rows):
        for j in range(columns):
            s = errorYList[j]
            tempList[j] = functionList[j][i]/s
        resultMatrix[i] = tempList
    return resultMatrix


# initialization of Matrix A
matrixA = matrixFunction(YList, errorYList)

# initialization and assignment of vector b
vectorB = np.empty(len(yMeasured))

for y, s, i in zip(yMeasured, errorYList, range(len(yMeasured))):
    vectorB[i] = y/s
    print y, s

# ((A transpose) dot (A))
tempA = np.dot(matrixA.T, matrixA)

# inverse of a matrix
tempB = tempA.I

# Dot Product of matrix and b vector
tempC = np.dot(matrixA.T, vectorB)

finalTemp = np.dot(tempC, tempB)
vectorA = np.array(finalTemp)[0]

for Yvals in YList:
    plt.plot(xList, Yvals)
