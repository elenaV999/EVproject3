#READ ALL FILES WITH PARTICLE POSITION DATA  
#->  writes file with cm position for the 2 galaxies at each timestep

import numpy as np
import matplotlib.pyplot as plt

N=int(2e4) #number of lines in the file
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



f = open('cm_positions','w')
f.write('#xcm1		ycm1		zcm1		xcm2		ycm2		zcm2		cm_separation\n')

m=504.027	#single particle mass
Mgal=m*1e4 #mass of 1 galaxy (Npart=1e4)



#generate file names
Nfiles=800
iOutInterval = 1 
n = np.arange(1, Nfiles+1, 1 ) #start, stop, step
nfile = np.char.zfill(n.astype('str'),3) 
for k in range(Nfiles):	#loop over all files
	fname=str(nfile[k])+".dat"  
	print(fname)
	x,y,z = readfast(fname, N)
	#x,y,z = np.genfromtxt(fname, dtype="float", comments="#", usecols=(0,1,2), unpack=True)
	
	x1, x2 = np.split(x,2)
	y1, y2 = np.split(y,2)
	z1, z2 = np.split(z,2)

	xcm1=sum(m*x1)/Mgal	
	ycm1=sum(m*y1)/Mgal
	zcm1=sum(m*z1)/Mgal 

	xcm2=sum(m*x2)/Mgal
	ycm2=sum(m*y2)/Mgal
	zcm2=sum(m*z2)/Mgal
	
	separation=np.sqrt((xcm1-xcm2)**2+(ycm1-ycm2)**2+(zcm1-zcm2)**2)
	
	f.write(str(xcm1)+'	'+str(ycm1)+'	'+str(zcm1)+'	'+str(xcm2)+'	'+str(ycm2)+'	'+str(zcm2)+'	'+str(separation)+'\n')

