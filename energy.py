#check energy conservation

import numpy as np
import matplotlib.pyplot as plt

fname="out.log"   #radial profiles of some quantities, projected on the plane orthogonal to the angular momentum
t, E, T, U, lx, ly, lz = np.genfromtxt(fname, dtype="float", comments="#", usecols=(0,2,3,4, 6,7,8), unpack=True)

dE=np.zeros(len(E))
for i in range(len(E)-1): 
	dE[i]=(E[i+1]-E[i])/E[i]
	
l = np.sqrt(lx**2+ly**2+lz**2)
dl=np.zeros(len(l))
for i in range(len(l)-1): 
	dl[i]=(l[i+1]-l[i])/l[i]

	

fig, ax = plt.subplots(2 ,2 , figsize=(7,6))
ax[0][0].plot(t , dE , marker='.',markersize='5',linestyle=':', color='b') 
ax[0][0].set_title('energy error', fontsize=13)
ax[0][0].set_ylabel('dE/E', fontsize=12)
ax[0][0].set_xlabel('t [Gyr]', fontsize=12)


ax[0][1].plot(t , dl , marker='.',markersize='5',linestyle=':', color='b') 
ax[0][1].set_title('momentum error', fontsize=13)
ax[0][1].set_ylabel('dL/L', fontsize=12)
ax[0][1].set_xlabel('t [Gyr]', fontsize=12)


ax[1][0].plot(t , T , marker='.',markersize='5',linestyle=':', color='b') 
ax[1][0].set_title('kinetic energy', fontsize=13)
ax[1][0].set_ylabel('T', fontsize=12)
ax[1][0].set_xlabel('t [Gyr]', fontsize=12)


ax[1][1].plot(t , U , marker='.',markersize='5',linestyle=':', color='b') 
ax[1][1].set_title('potental energy', fontsize=13)
ax[1][1].set_ylabel('U', fontsize=12)
ax[1][1].set_xlabel('t [Gyr]', fontsize=12)
fig.tight_layout() 
plt.show()











