import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import matplotlib

# Note: the plotting tools are a bit of a mess at the moment and are in dire need of
# an optimization and cleaning pass

def PlotFreesurfaces(U,m,H=400,radius=1.0,scalefactor=1.0):
    pl = [None]*len(U)	
    for i in range(0,len(U)):
        u = U[i]
        plt.figure(i)
        PlotFreesurface(i,radius,h=H,f=(lambda rad,theta: scalefactor*np.cos(m*theta)*u(rad,0)))
    plt.show(block=False)
    raw_input("Press enter to continue.")
    plt.close("all")

def PlotHeightProfiles(U,m,H=400,radius=1.0,scalefactor=1.0):
    figs = [0]*len(U)
    for i in range(0,len(U)):
        u = U[i]
        plt.figure(i)
        PlotProfile(radius,h=H,f=(lambda rad,theta: scalefactor*np.cos(m*theta)*u(rad,0)))
    plt.show(block=False)
    raw_input("Press enter to continue.")
    plt.close("all")

def DualPlots(U,m,H=400,radius=1.0,container='cylinder'):
    figs = [0]*len(U)  
    for i in range(0,len(U)):
	u = U[i]
	fig = plt.figure(i,figsize=(3,3))
	(sf,min,max) = ComputeScalingFactor(u,m,container)
	DualPlot(u,m,i,fig,scalingfactor=sf,umin=min,umax=max,domaintype=container,h=H)
    plt.show(block=False)
    raw_input("Press enter to continue.")
    plt.close("all")


def ComputeScalingFactor(u,m,container="cylinder",h=100):
    func = np.vectorize(u)
    U = np.append((-1.0)**m*func(np.linspace(1.0,0,h),np.linspace(0,0,h)),
		            func(np.linspace(0,1.0,h),np.linspace(0,0,h))) 
    umin = np.abs(np.amin(U))
    umax = np.abs(np.amax(U))
    umm = max(umin,umax)
    if container == "cylinder":
    	return (0.4/umm,umin,umax)
    if container == "halfcircle":
        return (0.4/umm,umin,umax)
    if container == "tallglass":
	return (0.5/umm,umin,umax)
    if container == "dentedglass":
	return (0.5/umm,umin,umax)
    

def DualPlot(u,m,n,fig,scalingfactor=1.0,umin=-1.0,umax=1.0,domaintype="cylinder",h=100,radius=1.0,lw=3):
    func = np.vectorize(u)
    U = np.append((-1.0)**m*func(np.linspace(1.0,0,h),np.linspace(0,0,h)),
		            func(np.linspace(0,1.0,h),np.linspace(0,0,h))) 
      
 

    if domaintype == "halfcircle":
	X = np.linspace(-1,1,h)
	f = lambda x: -1.0*np.sqrt(1.0-x**2)
	plt.plot(X,f(X),linewidth=lw,color='black',solid_capstyle='round')
	
	[A,B] = np.meshgrid(np.linspace(-1,1,200),np.linspace(0,-1,200))
        

    if domaintype == "cylinder":
	X = np.linspace(-1,1,h)
	plt.vlines([-1,1],0,[-1,-1],linewidth=lw,color='black')
	plt.hlines([-1],-1,[1],linewidth=lw,color='black')
	
	[A,B] = np.meshgrid(np.linspace(-1,1,200),np.linspace(0,-1,200))
        

    if domaintype == "tallglass":
	rad = np.sqrt(2.0)
	X = np.linspace(-1.0*rad,1.0*rad,h)
        f = lambda x: -1.0*np.sqrt(rad**2-x**2) - 1.0
	plt.plot(X,f(X),linewidth=lw,color='black')
	X = np.linspace(-1.0*rad,-1.0,20)
	f = lambda x: 1.0*np.sqrt(rad**2-x**2) - 1.0
	plt.plot(X,f(X),linewidth=lw,color='black')
 	X = np.linspace(1.0,1.0*rad,20)
	plt.plot(X,f(X),linewidth=lw,color='black')

	[A,B] = np.meshgrid(np.linspace(-rad,rad,200),np.linspace(0,-1-rad,200))

    if domaintype == "dentedglass":
	rad = np.sqrt(2.0)
	X = np.linspace(-1.0*rad,1.0*rad,200)
	f = lambda x: -1.0*np.sqrt(rad**2-x**2) - 1.0
	plt.plot(X,f(X),linewidth=lw,color='black')

	X = np.linspace((1.175875602-0.2772343382),1.306562965,20)
	f = lambda x: -1.0*np.sqrt(0.2772343382**2 - (x-1.175875602)**2) + 0.7856949583 - 1.0
	plt.plot(X,f(X),linewidth=lw,color='black')
	X = np.linspace((1.175875602-0.2772343382),1.0,20)
	f = lambda x: np.sqrt(0.2772343382**2 - (x-1.175875602)**2) + 0.7856949583 - 1.0
	plt.plot(X,f(X),linewidth=lw,color='black')

	X = np.linspace(-1.306562965,(-1.175875602+0.2772343382),20)
	f = lambda x: -1.0*np.sqrt(0.2772343382**2 - (x+1.175875602)**2) + 0.7856949583 - 1.0
	plt.plot(X,f(X),linewidth=lw,color='black')
	X = np.linspace(-1.0,(-1.175875602+0.2772343382),20)
	f = lambda x: np.sqrt(0.2772343382**2 - (x+1.175875602)**2) + 0.7856949583 - 1.0
	plt.plot(X,f(X),linewidth=lw,color='black')

	X = np.linspace(1.306562965,np.sqrt(2.0),20)
	f = lambda x: 1.0*np.sqrt(rad**2-x**2) - 1.0
	plt.plot(X,f(X),linewidth=lw,color='black')
	X = np.linspace(-np.sqrt(2.0),-1.306562965,20)
	plt.plot(X,f(X),linewidth=lw,color='black')

	[A,B] = np.meshgrid(np.linspace(-rad,rad,200),np.linspace(0,-1-rad,200))	
	
    

    T = np.append(np.linspace(-1.0,0,h),np.linspace(0,1.0,h))

       
    try:
        u.set_allow_extrapolation(True)
    except:
	pass

    
    u0 = -scalingfactor*u(1.0,0)
    u1 = scalingfactor*u(1.0,0)

    # Plot the end bars and zero lines
    

    mask = A < 0
    C = func(abs(A),B)
    C[mask] = (-1.0)**m * C[mask]

    if domaintype=="halfcircle":
        C[A**2+B**2 > 1] = np.NaN

	plt.contourf(A,B,C,100,cmap='RdBu',vmin=-max([umin,umax]),vmax=max([umin,umax]))
        plt.ylim([-1.25,1.25])
        plt.xlim([-1.25,1.25])
	plt.vlines([-1,1],-scalingfactor*max([umin,umax])+(5.0/8.0),[scalingfactor*max([umin,umax])+(5.0/8.0),scalingfactor*max([umin,umax])+(5.0/8.0)],color='black',linewidth=lw)
    	plt.hlines([0,(5.0/8.0)],[-1,-1],[1,1],linestyle='dotted')
    	plt.plot(T,scalingfactor*U+(5.0/8.0),linewidth=2,color='black')
	plt.xticks([-1.25,0,1.25])
	plt.yticks([-1.25,0,1.25])
	plt.tick_params(axis='both',which='major',labelsize=8)
	#plt.savefig('plots/bowl_m'+str(m)+'_n'+str(n)+'.png')

    if domaintype=="tallglass":
	C[A**2+(B+1)**2 > np.sqrt(2.0)**2] = np.NaN
	plt.contourf(A,B,C,100,cmap='RdBu',vmin=-max([umin,umax]),vmax=max([umin,umax]))
        plt.ylim([-2.5,1.5])
        plt.xlim([-2,2])

	X = np.linspace(-2**0.5,2**0.5,h) 
	D = func(1.0/(2**0.5) * abs(X), np.zeros(h))
	D[X<0] = (-1.0)**m * D[X<0]
	plt.vlines([-2**0.5,2**0.5],-scalingfactor*max([umin,umax])+0.75,[scalingfactor*max([umin,umax])+0.75,scalingfactor*max([umin,umax])+0.75],color='black',linewidth=lw)
    	plt.hlines([0,0.75],[-1,-2**0.5],[1,2**0.5],linestyle='dotted')
    	plt.plot(X,scalingfactor*D+0.75,linewidth=2,color='black')
	plt.xticks([-2,0,2])
	plt.yticks([-2.5,0,1.5])
	plt.tick_params(axis='both',which='major',labelsize=8)
	#plt.savefig('plots/tallgass_m'+str(m)+'_n'+str(n)+'.png')


    if domaintype=="dentedglass":
	c1 = [1.175875602,-0.2143050417]
	c2 = [-1.175875602,-0.2143050417]
	radius1 = 0.2772343382
	C[A**2+(B+1)**2 > np.sqrt(2.0)**2] = np.NaN
	C[(A-c1[0])**2+(B-c1[1])**2 < radius1**2] = np.NaN
	C[(A-c2[0])**2+(B-c2[1])**2 < radius1**2] = np.NaN
        plt.contourf(A,B,C,100,cmap='RdBu',vmin=-max([umin,umax]),vmax=max([umin,umax]))
        plt.ylim([-2.5,1.5])
        plt.xlim([-2,2])

	X = np.linspace(-2**0.5,2**0.5,h) 
	D = func(1.0/(2**0.5) * abs(X), np.zeros(h))
	D[X<0] = (-1.0)**m * D[X<0]
	plt.vlines([-2**0.5,2**0.5],-scalingfactor*max([umin,umax])+0.75,[scalingfactor*max([umin,umax])+0.75,scalingfactor*max([umin,umax])+0.75],color='black',linewidth=lw)
    	plt.hlines([0,0.75],[-1,-2**0.5],[1,2**0.5],linestyle='dotted')
	plt.plot(X,scalingfactor*D+0.75,linewidth=2,color='black')
	plt.xticks([-2,0,2])
	plt.yticks([-2.5,0,1.5])
	plt.tick_params(axis='both',which='major',labelsize=8)
	#plt.savefig('plots/dentedglass_m'+str(m)+'_n'+str(n)+'.png')


    if domaintype=="cylinder":
        plt.contourf(A,B,C,100,cmap='RdBu',vmin=-max([umin,umax]),vmax=max([umin,umax]))
        plt.ylim([-1.25,1.25])
        plt.xlim([-1.25,1.25])
	plt.vlines([-1,1],-scalingfactor*max([umin,umax])+(5.0/8.0),[scalingfactor*max([umin,umax])+(5.0/8.0),scalingfactor*max([umin,umax])+(5.0/8.0)],color='black',linewidth=lw)
    	plt.hlines([0,(5.0/8.0)],[-1,-1],[1,1],linestyle='dotted')
    	plt.plot(T,scalingfactor*U+(5.0/8.0),linewidth=2,color='black')
	plt.xticks([-1.25,0,1.25])
	plt.yticks([-1.25,0,1.25])
	plt.tick_params(axis='both',which='major',labelsize=8)
	#plt.savefig('plots/cyl_m'+str(m)+'_n'+str(n)+'.png')


def AnimateDualPlot(u,m,n,fig,container="cylinder",movielength=1.5,H=400):
    Fps = 90
    #Nframes = movielength*Fps
    Nframes = 180
    Tj = lambda j: j/(1.0*Nframes)
    u.set_allow_extrapolation(True)
    (sf,min,max) = ComputeScalingFactor(u,m)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=Fps,bitrate=1800)

    plt.clf()
    def Update(nframe):
        plt.clf()
        DualPlot((lambda r,z: np.cos(2*np.pi*Tj(1.0*nframe))*u(r,z)),m,n,fig,scalingfactor=sf,umin=min,umax=max,domaintype="cylinder",h=H)
	print(str(nframe)+'/'+str(Nframes))

    anim = animation.FuncAnimation(fig,Update,frames=Nframes,interval=11.11)
    
    if container == "cylinder":    
	anim.save('cyl_m'+str(m)+'_n'+str(n)+'_anim.mp4',writer)
    elif container == "halfsphere":    
	anim.save('bowl_m'+str(m)+'_n'+str(n)+'_anim.mp4',writer)
    elif container == "tallglass":    
	anim.save('tallglass_m'+str(m)+'_n'+str(n)+'_anim.mp4',writer)
    elif container == "dentedglass":    
	anim.save('dentedglass_m'+str(m)+'_n'+str(n)+'_anim.mp4',writer)
    


# Helper Functions
def PlotProfile(radius,theta=0,h=10,f=(lambda x,y: 0)):
    R = np.append(np.linspace(radius,0,h),np.linspace(0,radius,h))			# Initialize sample point arrays
    P = np.append(np.linspace(theta+np.pi,theta+np.pi,h),np.linspace(theta,theta,h))

    T = np.linspace(-1,1,2*h)								# Initialize parameter array

    func = np.vectorize(f)				# Vectorize function of interest and sample at the points (R,P)
    W = func(R,P)

    plt.plot(T,W)					# Plot arbitrary parameter T:[-1,1] against samples array W

def PlotFreesurface(I,radius,h=10,f=(lambda x,y: 0)):
    rhorange = np.linspace(0,radius,h)			# Initialize sample point arrays
    thetarange = np.linspace(0,2*np.pi,h)
    R, P = np.meshgrid(rhorange, thetarange)

    func = np.vectorize(f)				# Vectorize function of interest and sample at the points (R,P)
    Z = func(R,P)
    X, Y = R*np.cos(P), R*np.sin(P)			# Convert to cartesian
    
    fig = plt.figure(I)					# Initialize figure and axes then plot surface
    ax = fig.gca(projection="3d")
    return ax.plot_surface(X,Y,Z,cmap='autumn')
