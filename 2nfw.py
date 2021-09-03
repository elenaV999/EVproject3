####################################### by M Mapelli on October 20 2020##########################################################
import scipy.integrate as integ
import scipy.special as special
import numpy as np
import matplotlib.pyplot as plt

npart=int(1e4)	#number of particles
nbins = 20001	#number of bins for sampling of probability (generate equally spaced values for x,v)

msun=1.989e33 #in grams
km=1.0e5 #in cm
parsec=3.085678e18 #in cm
yr=3.1536e7  #in seconds
G_=6.674e-08	#G in cm^3/(s^2*grams)


#scale units and constants for Nbody units
scalength=1e3*parsec	#1kpc in cm
scalemass=2.25e5*msun    #2.25e5 Msun  	in grams		--> set this way to get tscale=1 Gyr	, velocity=kpc/Gyr=0.94km/s


#check units: 
print('length scale:		1kpc = ',scalength, ' cm' )
print('mass scale: 		2.25e5 Msun = ', scalemass, ' grams')
scaletime=np.sqrt(scalength**3/(G_*scalemass))	#1Gyr in seconds (set to have G=1)
print('-> time scale = ','{:e}'.format(scaletime/yr), 'yr = 1 Gyr')		#=10^9 yr = 1Gyr

vscale = 1e-5*scalength/scaletime #1 Gyr/kpc in km/s 
print('-> velocity scale = ', vscale, 'km/s = 1 Kpc/Gyr ') 
G=G_*scalemass*scaletime**2/scalength**3	#G=1 in Nbody units
print('G in Nbody units:		', G)


################# set parameters  ####################
#galaxy parameters: (optimized for the Milky Way)
c1=12.    #NFW concentration 
r200_1=220. #virial radius [kpc] 
c2=12.    
r200_2=220. 

#orbit parameters in Nbody units
k=2 #CHOOSE WHICH VALUE OF THA ARRAY TO USE:
vtan_kms = np.array([26.,57., 92.])#in km/s	<-----------------------MODIFY THIS!!!
vtan = vtan_kms/vscale	#tangential velocity in kpc/Gyr
print('tangential velocity ', vtan)
vrad=120./vscale #radial velocity [kpc/Gyr]
print('radial velocity = 120 km/s = ', vrad, 'kpc/Gyr')
vrel=np.sqrt(vtan**2+vrad**2) #total velocity [kpc/Gyr]
print('total velocity', vrel)

sep=780.  #[kpc]	initial separation
alpha = np.arctan(vtan/vrad)  
b = sep*np.sin(alpha) #impact parameter (separation along y)
print(len(alpha), len(b))
x_sep = sep*np.cos(alpha)#separation along x

print('alpha', alpha)
print('b', b, 'kpc')
print('x_sep', x_sep, ' kpc')
print('check Rsep=780kpc', np.sqrt(b**2+x_sep**2))

b_=[0., 0.4*r200_1, 0.8*r200_1]
print('previous b =', b_, 'kpc')

#cosmological parameters
zvir=0.0           #virialization redshift of NFW halo
omegam=0.3075        #omega_matter
omegalambda=0.691   #omega_lambda 
Hubble=6.774e01 / (3.085678e18 * 1.e01) #Hubble constant at z=0 in s
Hubblevir= Hubble * np.sqrt((omegam * ((1+zvir)**3.))+omegalambda) #Hubble constant at z=zvir  in s 

#other constants
pi= np.pi
rad_pi = np.sqrt(pi)
pi2 = pi*pi
G=(G_ * scalemass/(scalength**3)) # G in kpc^3 10^-8msun^-1 s^-2

rhocrit=((3. * Hubblevir**2/(8. * pi * G)))  #critical density

#############################NFW functions (see http://web.pd.astro.it/mapelli/images/thesis.ps.gz) ###############################
def sigmaNFW(x):  # rho(r) * M(r)  r-dependence inside the sigma(r) integral: rho(r)*M(r)/r^2
    f=1./x**3
    f=f/(1+x)**2
    f=f*(-x/(1.+x) + np.log(1.+x))
    return f

def fvel(v):   # P(v/sigma) = int Maxwellian       P(v/sigma)	cumulative of the maxwellian f(v,r) at radius r , v=v/sigma*sqrt(2)
    return ((-0.5 * np.exp(-v*v) * v)+(0.25 * np.sqrt(3.14159) * special.erf(v)))



def makeNFW (concentr, raggio200 ):
	####################INITIALIZE MASS AND RADII#########################
	#dependent NFW parameters
	deltac=(2.0e+02/3.0e0) * concentr * concentr * concentr/(np.log(1.0+concentr)- concentr/(1.0+concentr))
	rho200=rhocrit * deltac       #NFW rho200
	rscale=(raggio200/concentr)   #NFW scale radius

	raggiomin=(0.0 * rscale)
	M200=(2.0e2 * rhocrit * 4. * pi * (raggio200**3.)/3.)
	MTOTALE=M200 # NFW mass
	
	x200     = concentr
	xmax     = 1. * x200


	#Monte carlo sampling for positions and velocities 
	################INITIALIZE SAMPLING ARRAYS##########################  
	x_list     = np.zeros(nbins)
	Mx_list    = np.zeros(nbins)
	Sigma_list = np.zeros(nbins)
	vs_list    = np.zeros(nbins)
	fvs_list   = np.zeros(nbins)

	vs_max = 10.   # neglect/exclude particles with speed > 10 sigma

	x_list[0]=0.0	
	vs_list[0] = 0.0
	for i in range(1,nbins):
	    x_list[i]=x_list[i-1]+1.0
	    vs_list[i] = vs_list[i-1]+1.0
	    #print(x_list[i])


	#################sample probability of positions####################### 
	x_list[:] = xmax * x_list[:] / (float(nbins)-1.)	#equally spaced values for x
	#print("Mx_list")
	Mx_list[:] = -x_list[:]/(1.+x_list[:]) + np.log(1.0+x_list[:])	#calculate M(x) [cumulative distr for the mass] proportional to P(x) [cumulative for the position]


	###############calculate velocity dispersion############################  for each equally spaced point x	
	p=integ.quad(sigmaNFW, (0.01*x_list[1]), (1.01*xmax), epsabs=1.49e-08, epsrel=1.49e-08, limit=50000)   #integrate the r-dependence

	Sigma_list[0] =p[0]

	for i in range (1,nbins): 
	    p= integ.quad(sigmaNFW, x_list[i], 1.01*xmax,  epsabs=1.49e-08, epsrel=1.49e-08, limit=5000)	
	    Sigma_list[i] = p[0]
	Sigma_list[:] = Sigma_list[:] * x_list[:] *((1.+x_list[:])**2)	
	Sigma_list[:] = Sigma_list[:] * (4. * pi * G * rho200 * (rscale**2))  #multipy by the needed factors ouside the integral
	Sigma_list[:] = np.sqrt(Sigma_list[:]) * (scalength/km)	

	###########sample probability of velocities P(v,r)##########################
	fvel0 = fvel(0.0)

	vs_list[:] = vs_max * vs_list[:] / (float(nbins)-1.) 	#equally spaced values for v/sigma
	fvs_list[:] = (4./rad_pi) * ( fvel(vs_list/np.sqrt(2.)) - fvel0 )	#P(v/sigma) = P(v,r) [cumulative maxwellian]
	 
	# ensure normalization to 1
	Mtot_halo= Mx_list[nbins-1] 
	fvs_max = fvs_list[nbins-1] 

	Mx_list[:] = Mx_list[:] / Mtot_halo	#P(x) = P(r) [cumulative for the position]  	x=r/rs
	fvs_list[:] = fvs_list[:] / fvs_max	#P(v/sigma(r)) = P(v,r) [maxwellian f(v,r) fumulative]

  	############GENERATE PARTICLE POSITIONS AND VELOCITIES BY INTERPOLATION###############    
	raggio=0.0
	r_seeds = np.zeros(npart)
	v_seeds = np.zeros(npart)
	x_part = np.zeros(npart)
	sigma_part = np.zeros(npart)
	vs_part = np.zeros(npart)

	x = np.zeros(npart)
	y = np.zeros(npart)
	z = np.zeros(npart)
	velx = np.zeros(npart)
	vely = np.zeros(npart)
	velz = np.zeros(npart)


	r_seeds = np.random.random(npart)	#ransom sample P(x)
	v_seeds = np.random.random(npart)	#ransom sample P(v,r)

	# interpolations
	x_part= np.interp(r_seeds,Mx_list,x_list)    #arguments: y=P(x) samples, P(x) calculated at every x, equally spaced x 
										#find positions x NFW-distributed
	sigma_part = np.interp(x_part, x_list, Sigma_list) #interpolate Sigma 	
							#arguments: corretly distributed x values, equally spaced x, sigma calculated at every equally spaced x
							#find sigma for a NFW distribution
	vs_part = np.interp(v_seeds,fvs_list,vs_list) #interpolate and invert fvs_list to get particle velocities
							#arguments: y = P(v) samples, P(v) calculated at every v, equally spaced v
							#find velocities v maxwell-distributed

    
	#########  generate angular coord and conversions to cartesian #############
	#POSITIONS
	raggio = x_part[:] * rscale		#angular coord (sperical symm -> uniform theta, phi)
	fi       = np.random.random(npart) * 2. * pi
	ctheta   = 2. * np.random.random(npart) - 1.
	theta    = np.arccos(ctheta)

	x=(raggio * np.sin(theta) * np.cos(fi))	# go cartesian
	y=(raggio * np.sin(theta) * np.sin(fi))
	z=(raggio * np.cos(theta))


	#VELOCITIES
	velocity = vs_part[:] * sigma_part[:]	
		
	#angular coord (velocity directon is isotropic ->  uniform theta, phi)
	vfi     = np.random.random(npart) * 2. * pi
	cvtheta = 2.*np.random.random(npart) - 1.
	vtheta  = np.arccos(cvtheta)

	velx=velocity * np.sin(vtheta) * np.cos(vfi) # go cartesian
	vely=velocity * np.sin(vtheta) * np.sin(vfi)
	velz=velocity * np.cos(vtheta)
		
	s=0.017 *rscale *(float(npart)/1e5)**-0.23  #dehnen 2001 SOFTENING  (see slides)
	softening = np.zeros(npart) + s
	
	
	return (M200, softening, x,y,z,velx,vely,velz)
	

#################################################### MAIN ##################################À

M200_1,s1,x1,y1,z1,velx1,vely1,velz1 = makeNFW(c1,r200_1)
M200_2,s2,x2,y2,z2,velx2,vely2,velz2 = makeNFW(c2,r200_2)

mass1 = np.zeros(npart)+ (M200_1/float(npart)) #single particle mass
mass2 = np.zeros(npart)+ (M200_2/float(npart)) 


#plt.scatter(x1,velx1,s=0.1) #
#plt.scatter(x2,velx2,s=0.1)
#plt.xlabel("x [kpc]")
#plt.ylabel("v$_x$ [kpc Gyr$^{-1}$]")
#plt.show()


###########################place in different POSITIONS
x1 += x_sep[k]/2 
x2 -= x_sep[k]/2
y1 += b[k]/2 
y2 -= b[k]/2
###########################add relative velocity
velx1 -= vrel[k]/2 
velx2 += vrel[k]/2


mass=np.concatenate((mass1, mass2))
s=np.concatenate((s1, s2))
x=np.concatenate((x1, x2))
y=np.concatenate((y1, y2))
z=np.concatenate((z1, z2))
velx=np.concatenate((velx1, velx2))
vely=np.concatenate((vely1, vely2))
velz=np.concatenate((velz1, velz2))

#Ek=0.5*mass[0]*(velx**2+vely**2+velz**2).sum()
#print(Ek)


######################## set CM=0, vCM=0 #########################
xcm = sum(mass*x)/(M200_1+M200_2)
ycm = sum(mass*y)/(M200_1+M200_2)
zcm = sum(mass*z)/(M200_1+M200_2)
rcm = np.sqrt(xcm**2+ycm**2+zcm**2) 
#print('cm = ', xcm, ycm, zcm, '\nr_cm = ', rcm)
vxcm = sum(mass*velx)/(M200_1+M200_2)
vycm = sum(mass*vely)/(M200_1+M200_2)
vzcm = sum(mass*velz)/(M200_1+M200_2)
vcm = np.sqrt(vxcm**2+vycm**2+vzcm**2) 
#print('vcm = ', vxcm, vycm, vzcm, '\nv_cm = ', vcm)
x-=xcm
y-=ycm
z-=zcm
velx-=vxcm
vely-=vycm
velz-=vzcm

#Ek=0.5*mass[0]*(velx**2+vely**2+velz**2).sum()
#print(Ek)


######################### OUTPUT FILES  ###################################
#recap of some parameters
filename = '0_vtan'+str(vtan_kms[k])+'.txt'
par_file =open(filename, 'w') 
par_file.write('Npart =	'+str(npart)+'\n')
par_file.write('softening = '+str(s1[3])+'\n')
par_file.write('rhocrit = '+str(rhocrit)+'  [N.B.units]\n\n')
par_file.write('galaxy parameters:\nc =		'+str(c1)+'\n')
par_file.write('R200 =	'+str(r200_1)+' [kpc]\n')
par_file.write('M200 =	'+str(M200_1)+' [N.B.units]	'+str(M200_1*2.25e5)+' M_sun\n\n')
par_file.write('orbit parameters:\nv_tan = '+str(vtan_kms[k])+' [km/s]	'+str(vtan[k])+' [kpc/Gyr]\n'  )
par_file.write('v_rad = 120 [km/s]	'+str(vrad)+' [kpc/Gyr]\n'  )
par_file.write('v = '+str(np.sqrt(120**2+vtan_kms[k]**2))+' [km/s]	'+str(vrel[k])+' [kpc/Gyr]'  )
par_file.write('Rsep =	'+str(sep)+' [kpc] \nb =		'+str(b[k])+' [kpc]\nx_sep = '+str(x_sep[k])+'[kpc]\n\n')
par_file.close()


#write initial conditions in tipsy ascii format
f=open("NFW_000000.tip","w")
f.write(str(2*npart)+' '+str(0)+' '+str(0)+'\n')	#npart=total n° of particles, ngas=0, nstar=0
f.write(str(3)+'\n')		#ndimension = 3 
f.write(str(0.0)+'\n')		#time=0 (initial conditions: we are at initial time)
for i in range(2*npart):
    f.write(str(mass[i])+'\n')
print("particle mass (same for all)")
print("mass[0], mass[npart-1], npart", mass[0], mass[2*npart-1], npart)
for i in range(2*npart):			#particle positions		
    f.write(str(x[i])+'\n')
#print("x[0], x[npart-1],npart", x[0], x[npart-1],npart)
for i in range(2*npart):
    f.write(str(y[i])+'\n')
for i in range(2*npart):
    f.write(str(z[i])+'\n')
for i in range(2*npart):
    f.write(str(velx[i])+'\n')	#particle velocities
for i in range(2*npart):
    f.write(str(vely[i])+'\n')
for i in range(2*npart):
    f.write(str(velz[i])+'\n')
for i in range(2*npart):
    f.write(str(s[i])+'\n')  
for i in range(2*npart):
    f.write(str(0.0)+'\n')	#potential energy (set to zero for the initial conditions)
f.close()


############### PLOT #############		x_part	sigma_part



#x-y plane
plt.scatter(x,y,s=0.1)
plt.xlabel("x [kpc]")
plt.ylabel("y [kpc]")
plt.show()

'''
#vx vy
plt.scatter(vely,velx,s=0.1)
plt.xlabel("v$_x$ [kpc Gyr$^{-1}$]")
plt.ylabel("v$_y$ [kpc Gyr$^{-1}$]")
plt.show()

plt.scatter(x,velx,s=0.1)
plt.xlabel("x [kpc]")
plt.ylabel("v$_x$ [kpc Gyr$^{-1}$]")
plt.show()


#sigma
plt.scatter(x_part, sigma_part,s=0.1)
plt.xlabel("sigma ")
plt.ylabel("x_part")
plt.show()


#NFW density distribution
Np=100
r=np.sqrt(x**2+y**2+z**2)
binr=np.linspace(min(r),max(r),Np)
#binr=np.logspace(np.log10(min(r)),np.log10(max(r)),50)

binm=np.zeros([Np-1])
rho=np.zeros([Np-1])
NFW=np.zeros([Np-1])
binrm=np.zeros([Np-1])

for i in range(1,len(binr)):
    binrm[i-1]=0.5*(binr[i-1]+binr[i])

for i in range(0,len(r)):
    for j in range(1,len(binr)):
        if((r[i]>=binr[j-1]) and (r[i]<binr[j])):
            binm[j-1]+=mass[i]
for i in range(0,len(binr)-1):
    rho[i]=4.*np.pi/3.*(binr[i+1]**3-binr[i]**3)
    rho[i]=binm[i]/rho[i]
    NFW[i]=rhocrit*deltac/((binrm[i]/rscale)*(1+binrm[i]/rscale)**2)

    
plt.plot(binrm,rho)
plt.plot(binrm,NFW)
plt.legend(["ICs", "Theory"])
plt.xscale("log")
plt.yscale("log")
plt.xlabel("r [kpc]")
plt.ylabel("$\\rho{}$ [Nbody units]")

plt.show()
'''


