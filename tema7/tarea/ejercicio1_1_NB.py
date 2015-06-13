import Image
import numpy as np
import matplotlib.pyplot as plt

import sys

def get_cabeza(INPUT, down_i, up_i, down_j, up_j):
        J= Image.open(INPUT)
        a=np.asarray(J,dtype=np.float32)

        cabeza = a[down_i:up_i, down_j:up_j]
        output = Image.fromarray(cabeza.astype(np.uint8))

        return [a, output]

if __name__ == "__main__":
	if len(sys.argv) == 6:
		IMAGE  = sys.argv[1]
		down_i = int(float(sys.argv[2]))
		up_i   = int(float(sys.argv[3]))
		down_j = int(float(sys.argv[4]))
		up_j   = int(float(sys.argv[5]))
		OUTPUT = "output_cameraman.jpg"
		output = get_cabeza(IMAGE, down_i, up_i, down_j, up_j)
		print "Matrix: %s" %(str(output[0]))
		output[1].save(OUTPUT)
		print "Output saved in: %s" %(OUTPUT)
		
	else:
		OUTPUT = "output_cameraman.jpg"
		IMAGE  = "cameraman_512.tif" 
		down_i = 70
		up_i   = 170
		down_j = 180
		up_j   = 280

		output = get_cabeza(IMAGE, down_i, up_i, down_j, up_j)
		print "Matrix: %s" %(str(output[0]))
		output[1].save(OUTPUT)
		print "Output saved in: %s" %(OUTPUT)
