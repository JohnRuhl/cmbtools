import numpy as np

test1 = np.array([[range(10), "freqs"], [range(15), "freqs"]])
test2 = np.array([range(20), "freqs"])

def testfunc(inputs):
    freqs = [0, 1, 2]
    evaln = [0]*len(freqs)  # each freq is an elt
    ind = np.where(inputs == "freqs")[0][0]  # index of freqs  # *** FIX ***
    print("ind", ind)
    for i, freq in enumerate(freqs):  # frequencies
        inputs[ind] = freq  # substitute freq
        evaln[i] = inputs[ind]
    inputs[ind] = "freqs"  # resetting inputs
    return evaln, inputs

def testfunc2(inputs):
    freqs = [0, 1, 2]
    ind = np.where(inputs == "freqs")
    print("ind", ind)
    for i, freq in enumerate(freqs):  # frequencies
        inputs[ind] = freq  # substitute freq
        evaln[i] = inputs[ind]


test3 = np.array([[range(10), "freqs"], [range(15), "freqs"]])
ind = np.where(test2 == "freqs")
test2[ind] = 999
print("no indices", test2)

evaln, inputs = testfunc(test2)
print("outputs:\n", evaln, inputs)
