import numpy as np
Nfiles=800
iOutInterval = 1 
n = np.arange(0, Nfiles+1, 1 ) #start, stop, step
nfile = np.char.zfill(n.astype('str'),3) #add zeros in front: array, n° of digits

f=open("ma_pos","w")
f.write('ma_pos\nshell mkdir positions\n')
for k in nfile:
	f.write('openbinary out.000'+k+'\nloads 1\nzall\nxyzpoints all ./positions/'+k+'.dat\nclosebinary\n')
f.write('end\n')
f.close()	



