import Image
import numpy as np
import matplotlib.pyplot as plt

import sys

def get_centro():
	I=Image.open("lena.jpg")
	I1=I.convert('L') 
	a=np.asarray(I1,dtype=np.float32)	

	b = np.zeros(shape=(512,512))
	[c_x, c_y] = [512/2, 512/2]
	for i in range(0,512):
		    for j in range(0,512):
		        b[i][j] = (i-c_x)**2+(j-c_y)**2 < 150**2
	result = a*b

	return result

def get_centro_m():
	I=Image.open("lena.jpg")
	I1=I.convert('L') 
	a=np.asarray(I1,dtype=np.float32)	

	b = np.zeros(shape=(512,512))
	[c_x, c_y] = [512/2, 512/2]
	for i in range(0,512):
		    for j in range(0,512):
		        if (i-c_x)**2+(j-c_y)**2 < 150**2:
				b[i][j]=1
			else:
				b[i][j]=0.5
	result = a*b

	return result

if __name__ == "__main__":
	if len(sys.argv) == 1:
		OUTPUT = "output_lena.jpg"
		OUTPUT_M = "output_lena_half.jpg"

		output = get_centro()
		Image.fromarray(output.astype(np.uint8)).save(OUTPUT)
		print "Output saved in: %s" %(OUTPUT)

		output_m = get_centro_m()
		Image.fromarray(output_m.astype(np.uint8)).save(OUTPUT_M)
		print "Output saved in: %s" %(OUTPUT_M)
		
	else:
		print "Usage: ejercicio1_2_NB.py" 
