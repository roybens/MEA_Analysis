import numpy as np
import scipy
def unweighted_entropy_score(arr1: np.array, arr2: np.array):
    if len(arr1) != len(arr2):
        raise ValueError('Length of two arrays should be the same')
        return
    arr1 = normalize_arr(arr1)
    arr2 = normalize_arr(arr2)
    sum = 0
    for i in range(len(arr1)):
        sum += entropy_match(arr1[i], arr2[i])
    return sum/2
def weigthed_entropy_score(arr1: np.array, arr2: np.array):
    if len(arr1) != len(arr2):
        raise ValueError('Length of two arrays should be the same')
        return
    arr1 = (normalize_arr(weight_intensity(arr1)))
    arr2 = (normalize_arr(weight_intensity(arr2)))
    sum = 0
    for i in range(len(arr1)):
        sum += entropy_match(arr1[i], arr2[i])
    return sum/2
# below are helper functions
def entropy_func(x: float):
    return x * np.log2(x)
def weight_intensity(arr):
    S = scipy.stats.entropy(arr)
    # print(S)
    if S >= 3:
        return arr
    else:
        w = 0.25+S*0.25
        return np.power(arr, w)

def entropy_match(m: float,
                  n: float):
    return entropy_func(m+n) - entropy_func(m) - entropy_func(n)
def normalize_arr(arr: np.array):
    return arr/np.sum(arr)