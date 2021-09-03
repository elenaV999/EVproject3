import numpy as np
Nfiles=800
iOutInterval = 1 
n = np.arange(400, Nfiles+1, iOutInterval ) #start, stop, step
nfile = np.char.zfill(n.astype('str'),3) #add zeros in front: array, n° of digits

f=open("ma","w")
f.write('ma\nshell mkdir photos\n')
for k in nfile:
	f.write('openbinary out.000'+k+'\nloads 1\nzall\nhard movie ./photos/out.000'+k+'.xwd\nshell convert ./photos/out.000'+k+'.xwd ./photos/out.000'+k+'.png\nclosebinary\n')
f.write('end\n')
f.close()	


