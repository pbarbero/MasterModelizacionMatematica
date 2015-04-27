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
    for i in range(0,20):
        GradU=np.gradient(u)
        pupx=GradU[0] # partial x
        pupy=GradU[1] # partial y
        Grad_pupx=np.gradient(pupx) 
        Grad_pupy=np.gradient(pupy)
        p2up2x=Grad_pupx[0]
        p2up2y=Grad_pupy[1]
        if np.linalg.norm(GradU) == 0:
            lap=(p2up2x+p2up2y)/(sqrt(eps**2+norm(GradU)**2))
        else:
            lap=(p2up2x+p2up2y)/np.linalg.norm(GradU)
        u = u+tau*(lap+lbda*(u0-u))
    
    u1=255*(u-u.min()) / (u.max()-u.min())
    
    return u1

def yaroslavsky(b0, h, rho):
    b=np.copy(b0)
    for i in range(0,256):
        for j in range(0,256):
            iMin=max(i - rho, 1)
            iMax=min(i + rho, 256)
            jMin=max(j - rho, 1)
            jMax=min(j + rho, 256)
            I=b0[iMin:iMax,jMin:jMax]
            H=np.exp(-(I-b0[i][j])**2/h**2)
	    print sum(H*I)/sum(H)
            b[i][j]=sum(H*I)/sum(H)
    return b

def PSNR(u0, u):
    MSE = ((u0 - u) ** 2).mean(axis=None)
    # MSE = 0
    # for i in range(0,256):
    #     for j in range(0,256):
    #         MSE = MSE + (u0[i][j] - u[i][j])**2
    # MSE = MSE/256**2
    return 20*log10(u.max()/sqrt(MSE))

if __name__ == "__main__":
    	I = Image.open("boat_256.png")
    	a=np.asarray(I,dtype=np.float32)
	a0 = aplicar_ruido(a)
	a1 = np.copy(a)
	a2 = np.copy(a)

	OUTPUT_VT  = "output_boat_restaurada_VT.png"
	output_vt  = variacion_total(a1, 0.1, 0.75, 1*math.exp(-6))
	Image.fromarray(output_vt.astype(np.uint8)).save(OUTPUT_VT)	
	print "Output Boat Restaurada saved in: %s" %(OUTPUT_VT)

	OUTPUT_YA  = "output_boat_restaurada_YA.png"
	output_ya  = yaroslavsky(a2, 30, 30)
	Image.fromarray(output_ya.astype(np.uint8)).save(OUTPUT_YA)	
	print "Output Boat Restaurada saved in: %s" %(OUTPUT_YA)

	print "PSNR asociado a filtro Variacion Total: "
	print PSNR(a, output_vt)
	print PSNR(a, output_ya)


	
