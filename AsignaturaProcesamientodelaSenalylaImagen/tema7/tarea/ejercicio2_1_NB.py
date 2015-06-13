import Image
import numpy as np
import matplotlib.pyplot as plt
import math
from numpy import linalg as LA


def aplicar_ruido(uc):
    u0=uc+12.*np.random.standard_normal(uc.shape) # 0.04705882352 error
    u0=u0.clip(min=0, max=255)

    return u0

def variacion_total(u0, lbda, tau, eps):
    u=np.copy(u0)
    #u=u0
    for i in range(0,20):
        GradU=np.gradient(u)
        pupx=GradU[0]
        pupy=GradU[1]
        coefx = pupx / np.sqrt(eps**2 + pupx**2 + pupy**2)
        coefy = pupy / np.sqrt(eps**2 + pupx**2 + pupy**2)
        lap = np.gradient(coefx)[0] + np.gradient(coefy)[1]        
        u = u+tau*(lap+lbda*(u0-u))

    u=255*(u-u.min())/(u.max()-u.min())    
    return u

if __name__ == "__main__":
    	I = Image.open("boat_256.png")
    	uc = np.asarray(I,dtype=np.float32)
	u0 = aplicar_ruido(uc)
	u1 = u0

	OUTPUT  = "output_boat_restaurada_VT_ejer2.png"
	output  = variacion_total(u1, 0.1, 0.75, 1*math.exp(-6))
	Image.fromarray(output.astype(np.uint8)).save(OUTPUT)	
	print "Output Boat Restaurada saved in: %s" %(OUTPUT)
