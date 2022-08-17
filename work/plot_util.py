import numpy as np
import matplotlib.pyplot as plt    



def plot_line_data(filename,lw):

    dat = np.loadtxt(filename)

    x = dat[:,0]
    y = dat[:,1]

    plt.plot(x,y,'k',linewidth = lw)

    

def plot_2D_point_data(filename):
    

    dat = np.loadtxt(filename)

    x = dat[:,0]
    y = dat[:,1]

    f  = dat[:,2]
    fe = dat[:,3]

    plt.scatter(x,y,c=f,s = 25,cmap = 'RdBu')
#    plt.colorbar()
    plt.axis('equal')




def plot_2D_data(filename):
    

    mfile = np.loadtxt(filename)
    nx = int(mfile[0,0])
    ny = int(mfile[0,1])

    k = 0
    x = np.zeros((nx,ny))
    y = np.zeros((nx,ny))
    fun = np.zeros((nx,ny))
    
    for ix in range(nx):
        for iy in range(ny):
            k = k+1
            x[ix][iy]= mfile[k,0]
            y[ix][iy]= mfile[k,1]
            fun[ix][iy]= mfile[k,2]        



#    mm = np.max(abs(fun))
#    mm = np.max(-(fun))


    plt.figure()
    plt.pcolormesh(x,y,fun,shading='gouraud',cmap = 'RdBu')
    plt.colorbar()
#    plt.clim(-mm,mm)
    plt.axis('equal')



def plot_2D_vdata(filename):
    

    mfile = np.loadtxt(filename)
    nx = int(mfile[0,0])
    ny = int(mfile[0,1])

    k = 0
    x = np.zeros((nx,ny))
    y = np.zeros((nx,ny))
    fx = np.zeros((nx,ny))
    fy = np.zeros((nx,ny))
    fa = np.zeros((nx,ny))
    
    for ix in range(nx):
        for iy in range(ny):
            k = k+1
            x[ix][iy]= mfile[k,0]
            y[ix][iy]= mfile[k,1]
            fx[ix][iy]= mfile[k,2]
            fy[ix][iy]= mfile[k,3]
            fa[ix][iy]= np.sqrt(mfile[k,2]*mfile[k,2] + mfile[k,3]*mfile[k,3])



    mm = 10*np.max(fa)
            
    plt.figure()
    plt.pcolormesh(x,y,fa,shading='gouraud',cmap = 'jet')
    plt.colorbar()
#    plt.quiver(x,y,fx,fy)
    plt.quiver(x,y,fx,fy,units='xy',angles='xy',scale=mm,
               width = 0.004)

    plt.axis('equal')



def plot_2D_mesh_def(filename,ls,us):
    
    mfile = np.loadtxt(filename)
    nx = int(mfile[0,0])
    ny = int(mfile[0,1])

    k = 0
    x = np.zeros((nx,ny))
    y = np.zeros((nx,ny))
    ux = np.zeros((nx,ny))
    uy = np.zeros((nx,ny))
    ua = np.zeros((nx,ny))
    
    for ix in range(nx):
        for iy in range(ny):
            k = k+1
            x[ix][iy]= mfile[k,0]
            y[ix][iy]= mfile[k,1]
            ux[ix][iy]= mfile[k,2]
            uy[ix][iy]= mfile[k,3]
            ua[ix][iy]= np.sqrt(mfile[k,2]*mfile[k,2] + mfile[k,3]*mfile[k,3])



#    plt.figure()
#    plt.pcolormesh(x,y,ua,shading='gouraud',cmap = 'jet')

    for i in range(0,ny-1,ls):
        plt.plot(x[:,i]+us*ux[:,i],y[:,i]+us*uy[:,i],'k')

    for i in range(0,nx-1,ls):
        plt.plot(x[i,:]+us*ux[i,:],y[i,:]+us*uy[i,:],'k')

    plt.plot(x[:,0]+us*ux[:,0],y[:,0]+us*uy[:,0],'k',linewidth=3)
    plt.plot(x[:,ny-1]+us*ux[:,ny-1],y[:,ny-1]+us*uy[:,ny-1],'k',linewidth=3)

    plt.plot(x[0,:]+us*ux[0,:],y[0,:]+us*uy[0,:],'k',linewidth=3)
    plt.plot(x[nx-1,:]+us*ux[nx-1,:],y[nx-1,:]+us*uy[nx-1,:],'k',linewidth=3)

    plt.axis('equal')




def plot_2D_wavefront(filename):
    mfile = np.loadtxt(filename)
    nx = int(mfile[0,0])

    k = 0
    x = np.zeros((nx))
    y = np.zeros((nx))
    
    for ix in range(nx):
        k = k+1
        x[ix]= mfile[k,0]
        y[ix]= mfile[k,1]


    plt.plot(x,y,'.k',markersize=2)




def plot_2D_ray(filename):
    mfile = np.loadtxt(filename)
    nx = int(mfile[0,0])

    k = 0
    x = np.zeros((nx))
    y = np.zeros((nx))
    
    for ix in range(nx):
        k = k+1
        x[ix]= mfile[k,0]
        y[ix]= mfile[k,1]


    plt.plot(x,y,'k',linewidth=2.0)    


