import numpy as np

def array_center(arr):
	return np.array([(arr[i]+arr[i+1])/2 for i in range(len(arr[:-1]))])
