from __future__ import division
import Image
import numpy as np
import matplotlib.pyplot as plt

import sys

def degradado():
	I=Image.open("lena.jpg")
	I1=I.convert('L') 
	a=np.asarray(I1,dtype=np.float32)	
	x=np.linspace(0,1,512)
	mask = np.tile(x, (512, 1)).transpose()
	result = a*mask

	return result

if __name__ == "__main__":
	if len(sys.argv) == 1:
		OUTPUT = "output_lena_mask.jpg"

		output = degradado()
		Image.fromarray(output.astype(np.uint8)).save(OUTPUT)
		print "Output saved in: %s" %(OUTPUT)
	else:
		print "Usage: ejercicio1_2_NB.py" 
