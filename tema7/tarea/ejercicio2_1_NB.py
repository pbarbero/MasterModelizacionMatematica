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
    	uc = np.asarray(I,dtype=np.float32)
	u0 = aplicar_ruido(uc)
	u1 = u0

	OUTPUT  = "output_boat_restaurada.png"
	output  = restaurar(u1, 0.1, 0.75, 1*math.exp(-6))
	Image.fromarray(output.astype(np.uint8)).save(OUTPUT)	
	print "Output Boat Restaurada saved in: %s" %(OUTPUT)
