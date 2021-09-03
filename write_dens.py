#READ ALL FILES WITH PARTICLE POSITION DATA  
# writes file with positions of the density center of each galaxy, at every time-step

import numpy as np
import matplotlib.pyplot as plt
from statistics import mode



def readfast(fname, N):		#filename, number of lines in the file
    f=open(fname,"r")		
    (x,y,z) = (
        np.zeros(N,dtype="float"),	
        np.zeros(N,dtype="float"),
        np.zeros(N,dtype="float"))
    i=0
    for linetext in f:			 #loop to read all the lines in the file (linetext is a genetric line in the file)
        word_list = linetext.split()	#splits a line in elements
        x[i]=np.float(word_list[0])		 #fill numpy array with float data from the selected columns
        y[i]=np.float(word_list[1])
        z[i]=np.float(word_list[2])
        i=i+1	 #fill the next element of the arrays when we repeat the cycle
    f.close()
    return (x,y,z)


Npart=2e4
N=int(Npart)

Np= 20  #bin number
box = 100.  # histogram box coord (half side)
width=box/Np
print('bin n° = ', Np, 'box width = ', 2*box, ' kpc')
print('bin lenght = ', width, ' kpc') 


m=504.027	#single particle mass
Mgal=m*Npart/2 #mass of 1 galaxy (Npart=1e4)
Mtot=m*Npart   #total mass of the sysem


f = open('1dens_','w')
f.write('# time[Gyr]	x1[kpc]	y1	z1		x2	y2	z2 	\n')
f.write('# box half side = '+str(box)+' kpc\n# bin width = '+str(width)+' kpc\n')



start = 700
stop = 801  #will stop at the previous file
Nfiles=stop-start
n = np.arange(start, stop, 1 ) #start, stop, step
nfile = np.char.zfill(n.astype('str'),3) 
for k in range(Nfiles):	#loop over all files
	fname=str(nfile[k])+".dat"  
	print(fname)
	x,y,z = readfast(fname, N)

	#set the origin in the cm
	xcm = sum(m*x)/Mtot
	ycm = sum(m*y)/Mtot
	zcm = sum(m*z)/Mtot
	rcm = np.sqrt(xcm**2+ycm**2+zcm**2) 
	x-=xcm
	y-=ycm
	z-=zcm

	x1, x2 = np.split(x,2) #divide particles of the two galaxies
	y1, y2 = np.split(y,2)
	z1, z2 = np.split(z,2)
	

	############################### HISTOGRAMS to find the density center #############################À
	hx1, binx1 = np.histogram(x1, 2*Np, range=(-box,box))  #hist data, bin edges
	hx2, binx2 = np.histogram (x2, 2*Np, range=(-box,box))
	ix1 = np.argmax(hx1)
	ix2 = np.argmax(hx2)

	hy1, biny1 = np.histogram(y1, 2*Np, range=(-box,box))
	hy2, biny2 = np.histogram (y2, 2*Np, range=(-box,box))
	iy1 = np.argmax(hy1)
	iy2 = np.argmax(hy2)

	hz1, binz1 = np.histogram(z1, 2*Np, range=(-box,box))
	hz2, binz2 = np.histogram (z2, 2*Np, range=(-box,box))
	iz1 = np.argmax(hz1)
	iz2 = np.argmax(hz2)

	posx1 = binx1[ix1]+0.5*width
	posx2 = binx2[ix2]+0.5*width
	posy1 = biny1[iy1]+0.5*width
	posy2 = biny2[iy2]+0.5*width
	posz1 = binz1[iz1]+0.5*width
	posz2 = binz2[iz2]+0.5*width

	f.write(str(n[k]/100)+'   '+str(posx1)+'   '+str(posy1)+'   '+str(posz1)+'   '+str(posx2)+'   '+str(posy2)+'   '+str(posz2)+'\n')

f.close()
