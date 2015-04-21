import Image
import numpy as np
import matplotlib.pyplot as plt
import math
from numpy import linalg as LA


def aplicar_ruido(uc):
    u0=uc+12.*np.random.standard_normal(uc.shape) # 0.04705882352 error
    u0=u0.clip(min=0, max=255)

    return u0

def restaurar(u, lbda, tau, eps):
    for i in range(0,20):
        GradU=np.gradient(u)
        pupx=GradU[0] # partial x
        pupy=GradU[1] # partial y
        Grad_pupx=np.gradient(pupx) 
        Grad_pupy=np.gradient(pupy)
        p2up2x=Grad_pupx[0]
        p2up2y=Grad_pupy[1]
        if LA.norm(GradU) == 0:
            lap=(p2up2x+p2up2y)/(sqrt(eps**2+LA.norm(GradU)**2))
        else:
            lap=(p2up2x+p2up2y)/LA.norm(GradU)
        u = u+tau*(lap+lbda*(u0-u))
    
    u1=255*(u-u.min()) / (u.max()-u.min())
    
    return u1

if __name__ == "__main__":
    	I = Image.open("boat_256.png")
    	uc=np.asarray(I,dtype=np.float32)
	u0 = aplicar_ruido(uc)
	u2 = u1 = u0

	OUTPUT_1 = "output_boat_restaurada1.jpg"
	output1  = restaurar(u1, 10, 0.01, 1*math.exp(-6))
	Image.fromarray(output1.astype(np.uint8)).save(OUTPUT_1)	
	print "Output Restaurada 1 saved in: %s" %(OUTPUT_1)

	OUTPUT_2 = "output_boat_restaurada2.jpg"
	output2  = restaurar(u2, 0.1, 0.75, 1*math.exp(-6))
	Image.fromarray(output2.astype(np.uint8)).save(OUTPUT_2)	
	print "Output Restaurada 2 saved in: %s" %(OUTPUT_2)
