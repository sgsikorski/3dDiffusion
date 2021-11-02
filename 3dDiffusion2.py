"""
3dDiffusion.py

Written by Scott Sikorski
Last revised Sept 29 2021

3d diffusion equation solver
Density slices are outputted to screen

"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plot
from mpl_toolkits import mplot3d
import time

def mapPos(v,dv,vN,adjust):
    return clip(0,int((v-adjust)/dv),vN-1)

def clip(minVal,val,maxVal):
    return(max(min(val,maxVal),minVal))

def linearWeighting(xp,ximo,xi,xipo):
    if   (xp>ximo and xp<=xi):
        return (xp-ximo)/(xi-ximo)
    elif (xp>xi and xp<=xipo):
        return (xipo-xp)/(xipo-xi)
    else:
        return 0.0

def mapPosIndex():
    return

def main():
    xN = 64
    yN = 64
    zN = 64
    N = xN * yN * zN
    ppcx = ppcy = ppcz = 4
    # If these change, need to adjust the distrubution for the pos
    dsmcppc1 = 27
    dsmcppc2 = 125
    dsmcppc3 = 1000
    deg = 3
    gamma = 5.0/3.0

    dt = 0.001
    tf = 0.010
    nsteps = int(tf/dt)

    xmin = ymin = zmin = adjust = -1.0
    xmax = ymax = zmax = 1.0
    r0 = 0.3
    
    xx = np.linspace(xmin,xmax,xN)
    yy = np.linspace(ymin,ymax,yN)
    zz = np.linspace(zmin,zmax,zN)
    dx = (xmax - xmin) / xN
    dy = (ymax - ymin) / yN
    dz = (zmax - zmin) / zN
    V = dx * dy * dz

    X1, Y1 = np.meshgrid(xx,yy)
    X2, Y2 = np.meshgrid(xx,zz)
    X3, Y3 = np.meshgrid(yy,zz)

    rho0 = 1.0
    mass1 = rho0*V/dsmcppc1
    mass2 = rho0*V/dsmcppc2
    mass3 = rho0*V/dsmcppc3
    nParticles1 = nParticles2 = nParticles3 = 0
    particleCounter = np.zeros(((xN,yN,zN))).astype(int)

    dens0Grid = np.zeros(((xN,yN,zN)))
    dens = np.zeros(((xN,yN,zN)))
    dens2 = np.zeros(((xN, yN, zN)))
    dens3 = np.zeros(((xN, yN, zN)))

    densGrid = np.zeros(((xN,yN,zN)))
    massGrid = np.zeros(((xN,yN,zN)))

    const = (2*1*dt)**0.5
    next =next2 = next3 = 0

    for i in range(xN):
        for j in range(yN):
            for k in range(zN):
                r=(xx[i]**2+yy[j]**2+zz[k]**2)**0.5
                if r < r0: densGrid[i,j,k]=1.0


    # fig = plot.figure()
    # ax = plot.axes(projection='3d')

    # xaxis, yaxis = np.meshgrid(xx, yy)
    # zaxis = np.ones_like(xaxis)
    # test = densGrid[:, :, 12]

    # norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
    # ax.contour3D(xaxis,yaxis,zaxis)
    # #ax.plot_surface(xaxis, yaxis, zaxis, facecolors=plot.cm.jet(norm(test)))
    # # ax.plot_surface(xaxis, yaxis, zaxis, facecolors=plot.cm.jet(norm(densGrid)))

    # m = plot.cm.ScalarMappable(cmap=plot.cm.jet, norm=norm)
    # m.set_array([])
    # plot.colorbar(m)
    # plot.show()




    for i in range(xN):
        for j in range(yN):
            for k in range(zN):
                r=(xx[i]**2+yy[j]**2+zz[k]**2)**0.5
                if r < r0:
                    dens0Grid[i,j,k] = 1.0
                    nParticles1 += dsmcppc1
                    nParticles2 += dsmcppc2
                    nParticles3 += dsmcppc3
    
    posx = np.zeros(nParticles1)
    posy = np.zeros(nParticles1)
    posz = np.zeros(nParticles1)
    posx2 = np.zeros(nParticles2)
    posy2 = np.zeros(nParticles2)
    posz2 = np.zeros(nParticles2)
    posx3 = np.zeros(nParticles3)
    posy3 = np.zeros(nParticles3)
    posz3 = np.zeros(nParticles3)
    
    dsmcStart = time.perf_counter()
    for i in range(xN):
        for j in range(yN):
            for k in range(zN):
                r=(xx[i]**2+yy[j]**2+zz[k]**2)**0.5
                if r < r0:
                    # Distribution among cell
                    for ii in range(3):
                        for jj in range(3):
                            for kk in range(3):
                                posx[next] = xx[i] + dx*ii/3 + dx/6
                                posy[next] = yy[j] + dy*jj/3 + dy/6
                                posz[next] = zz[k] + dz*kk/3 + dz/6
                                next+=1
                    for ii in range(5):
                        for jj in range(5):
                            for kk in range(5):
                                posx2[next2] = xx[i] + dx*ii/5 + dx/10
                                posy2[next2] = yy[j] + dy*jj/5 + dy/10
                                posz2[next2] = zz[k] + dz*kk/5 + dz/10
                                next2+=1
                    for ii in range(10):
                        for jj in range(10):
                            for kk in range(10):
                                posx3[next3] = xx[i] + dx*ii/10 + dx/20
                                posy3[next3] = yy[j] + dy*jj/10 + dy/20
                                posz3[next3] = zz[k] + dz*kk/10 + dz/20
                                next3+=1
                
    for it in range(nsteps):
        print(f'DSMC: {it}/{nsteps}')
        particleCounter = np.zeros(((xN,yN,zN))).astype(int)
        particleCounter2 = np.zeros(((xN,yN,zN))).astype(int)
        particleCounter3 = np.zeros(((xN,yN,zN))).astype(int)
        ran = np.random.normal(0, 1, nParticles1*3)
        ran2 = np.random.normal(0,1,nParticles2*3)
        ran3 = np.random.normal(0,1,nParticles3*3)
        for ip in range(nParticles1):
            posx[ip] += const*ran[ip*3]
            posy[ip] += const*ran[ip*3+1]
            posz[ip] += const*ran[ip*3+2]
            xNew = mapPos(posx[ip], dx, xN, adjust)
            yNew = mapPos(posy[ip], dy, yN, adjust)
            zNew = mapPos(posz[ip], dz, zN, adjust)
            particleCounter[xNew, yNew, zNew] += 1
        for ip in range(nParticles2):
            posx2[ip] += const*ran2[ip*3]
            posy2[ip] += const*ran2[ip*3+1]
            posz2[ip] += const*ran2[ip*3+2]
            xNew = mapPos(posx2[ip], dx, xN, adjust)
            yNew = mapPos(posy2[ip], dy, yN, adjust)
            zNew = mapPos(posz2[ip], dz, zN, adjust)
            particleCounter2[xNew, yNew, zNew] += 1
        for ip in range(nParticles3):
            posx3[ip] += const*ran3[ip*3]
            posy3[ip] += const*ran3[ip*3+1]
            posz3[ip] += const*ran3[ip*3+2]
            xNew = mapPos(posx3[ip], dx, xN, adjust)
            yNew = mapPos(posy3[ip], dy, yN, adjust)
            zNew = mapPos(posz3[ip], dz, zN, adjust)
            particleCounter3[xNew, yNew, zNew] += 1
        dens = particleCounter*mass1/V
        dens2 = particleCounter2*mass2/V
        dens3 = particleCounter3*mass3/V

    print(f'Total time for dsmc: {0}', time.perf_counter() - dsmcStart)
    xabs,xweight = np.polynomial.hermite.hermgauss(ppcx)
    yabs,yweight = np.polynomial.hermite.hermgauss(ppcy)
    zabs,zweight = np.polynomial.hermite.hermgauss(ppcz)
    const = (2*dt)**0.5

    qdsmcStart = time.perf_counter()
    for it in range(nsteps):
        if it%1==0:
            print(f'Step: {it}/{nsteps}')

        massGrid = np.zeros(((xN,yN,zN)))
        for icx in range(xN):
            for icy in range(yN):
                for icz in range(zN):

                    if densGrid[icx,icy,icz] != 0:
                        for ipx in range(ppcx):
                            xNew = xx[icx] + const * xabs[ipx]
                            for ipy in range(ppcy):
                                yNew = yy[icy] + const * yabs[ipy]
                                for ipz in range(ppcz):
                                    zNew = zz[icz] + const * zabs[ipz]

                                    massPart = densGrid[icx,icy,icz]*V*xweight[ipx]*yweight[ipy]*zweight[ipz]/sum(xweight)/sum(yweight)/sum(zweight)

                                    xCell = mapPos(xNew,dx,xN,adjust)
                                    yCell = mapPos(yNew,dy,yN,adjust)
                                    zCell = mapPos(zNew,dz,zN,adjust)

                                    if xCell>0 and xCell<xN-2 and yCell>0 and yCell<yN-2 and zCell>0 and zCell<zN-2:
                                        wx1 = linearWeighting(xNew,xx[xCell-2],xx[xCell-1],xx[xCell])
                                        wx2 = linearWeighting(xNew,xx[xCell-1],xx[xCell],xx[xCell+1])
                                        wx3 = linearWeighting(xNew,xx[xCell],xx[xCell+1],xx[xCell+2])

                                        wy1 = linearWeighting(yNew,yy[yCell-2],yy[yCell-1],yy[yCell])
                                        wy2 = linearWeighting(yNew,yy[yCell-1],yy[yCell],yy[yCell+1])
                                        wy3 = linearWeighting(yNew,yy[yCell],yy[yCell+1],yy[yCell+2])

                                        wz1 = linearWeighting(zNew,zz[zCell-2],zz[zCell-1],zz[zCell])
                                        wz2 = linearWeighting(zNew,zz[zCell-1],zz[zCell],zz[zCell+1])
                                        wz3 = linearWeighting(zNew,zz[zCell],zz[zCell+1],zz[zCell+2])
                                        if  (wz1!=0):
                                            if (wy1!=0):
                                                massGrid[xCell-1,yCell-1,zCell-1]+=massPart*wx1*wy1*wz1
                                                massGrid[xCell,yCell-1,zCell-1]+=massPart*wx2*wy1*wz1
                                                massGrid[xCell+1,yCell-1,zCell-1]+=massPart*wx3*wy1*wz1
                                            if (wy2!=0):
                                                massGrid[xCell-1,yCell,zCell-1]+=massPart*wx1*wy2*wz1
                                                massGrid[xCell,yCell,zCell-1]+=massPart*wx2*wy2*wz1
                                                massGrid[xCell+1,yCell,zCell-1]+=massPart*wx3*wy2*wz1
                                            if (wy3!=0):
                                                massGrid[xCell-1,yCell+1,zCell-1]+=massPart*wx1*wy3*wz1
                                                massGrid[xCell,yCell+1,zCell-1]+=massPart*wx2*wy3*wz1
                                                massGrid[xCell+1,yCell+1,zCell-1]+=massPart*wx3*wy3*wz1
                                        ##################################
                                        if (wz2!=0):
                                            if (wy1!=0):
                                                massGrid[xCell-1,yCell-1,zCell]+=massPart*wx1*wy1*wz2
                                                massGrid[xCell,yCell-1,zCell]+=massPart*wx2*wy1*wz2
                                                massGrid[xCell+1,yCell-1,zCell]+=massPart*wx3*wy1*wz2
                                            if (wy2!=0):
                                                massGrid[xCell-1,yCell,zCell]+=massPart*wx1*wy2*wz2
                                                massGrid[xCell,yCell,zCell]+=massPart*wx2*wy2*wz2
                                                massGrid[xCell+1,yCell,zCell]+=massPart*wx3*wy2*wz2
                                            if (wy3!=0):
                                                massGrid[xCell-1,yCell+1,zCell]+=massPart*wx1*wy3*wz2
                                                massGrid[xCell,yCell+1,zCell]+=massPart*wx2*wy3*wz2
                                                massGrid[xCell+1,yCell+1,zCell]+=massPart*wx3*wy3*wz2
                                        ###################################
                                        if (wz3!=0):
                                            if (wy1!=0):
                                                massGrid[xCell-1,yCell-1,zCell+1]+=massPart*wx1*wy1*wz3
                                                massGrid[xCell,yCell-1,zCell+1]+=massPart*wx2*wy1*wz3
                                                massGrid[xCell+1,yCell-1,zCell+1]+=massPart*wx3*wy1*wz3

                                            if (wy2!=0):
                                                massGrid[xCell-1,yCell,zCell+1]+=massPart*wx1*wy2*wz3
                                                massGrid[xCell,yCell,zCell+1]+=massPart*wx2*wy2*wz3
                                                massGrid[xCell+1,yCell,zCell+1]+=massPart*wx3*wy2*wz3
                                            if (wy3!=0):
                                                massGrid[xCell-1,yCell+1,zCell+1]+=massPart*wx1*wy3*wz3
                                                massGrid[xCell,yCell+1,zCell+1]+=massPart*wx2*wy3*wz3
                                                massGrid[xCell+1,yCell+1,zCell+1]+=massPart*wx3*wy3*wz3

        densGrid = massGrid / V
    
    print(f'Total time for QDSMC: {0}', time.perf_counter() - qdsmcStart)

    # Adjust hard-coded numbers for slices
    lineoutx = densGrid[:, int(yN/2), int(zN/2)]
    dlineoutx = dens[:, int(yN/2), int(zN/2)]
    dlineoutx2 = dens2[:, int(yN/2), int(zN/2)]
    dlineoutx3 = dens3[:, int(yN/2), int(zN/2)]

    fig, ax = plot.subplots(4, 3, figsize=(10,10))
    cf4 = ax[0][0].contourf(X1, Y1, 1.2*dens[:,:,int(zN/2)], 20, cmap='RdGy')
    cf5 = ax[0][1].contourf(X2, Y2, 1.2*dens[:,int(yN/2),:], 20, cmap='RdGy')
    cf6 = ax[0][2].contourf(X3, Y3, 1.2*dens[int(xN/2),:,:], 20, cmap='RdGy')
    ax[0][0].set_title(f'DSMC (ppc=27)  Density(x,y) @ z=0')
    ax[0][1].set_title(f'Density(x,z) @ y=0')
    ax[0][2].set_title(f'Density(y,z) @ x=0')
    cf4 = ax[1][0].contourf(X1, Y1, 1.2*dens2[:,:,int(zN/2)], 20, cmap='RdGy')
    cf5 = ax[1][1].contourf(X2, Y2, 1.2*dens2[:,int(yN/2),:], 20, cmap='RdGy')
    cf6 = ax[1][2].contourf(X3, Y3, 1.2*dens2[int(xN/2),:,:], 20, cmap='RdGy')
    ax[1][0].set_title(f'DSMC (ppc=125)')
    ax[1][1].set_title(f'DSMC (ppc=125)')
    ax[1][2].set_title(f'DSMC (ppc=125)')
    cf4 = ax[2][0].contourf(X1, Y1, 1.2*dens3[:,:,int(zN/2)], 20, cmap='RdGy')
    cf5 = ax[2][1].contourf(X2, Y2, 1.2*dens3[:,int(yN/2),:], 20, cmap='RdGy')
    cf6 = ax[2][2].contourf(X3, Y3, 1.2*dens3[int(xN/2),:,:], 20, cmap='RdGy')
    ax[2][0].set_title(f'DSMC (ppc=1000)')
    ax[2][1].set_title(f'DSMC (ppc=1000)')
    ax[2][2].set_title(f'DSMC (ppc=1000)')
    cf1 = ax[3][0].contourf(X1, Y1, densGrid[:,:,int(zN/2)], 20, cmap='RdGy')
    cf2 = ax[3][1].contourf(X2, Y2, densGrid[:,int(yN/2),:], 20, cmap='RdGy')
    cf3 = ax[3][2].contourf(X3, Y3, densGrid[int(xN/2),:,:], 20, cmap='RdGy')
    ax[3][0].set_title(f"QDS (ppc=4)")
    ax[3][1].set_title(f"QDS (ppc=4)")
    ax[3][2].set_title(f"QDS (ppc=4)")
    cax = plot.axes([0.925,0.1,0.025,0.7])
    fig.colorbar(cf1, cax=cax)
    plot.show()

    fig, ax = plot.subplots(1, 1, figsize=(10,10))
    ax.plot(xx,lineoutx, color='black',label='QDSMC ppc=4')
    ax.plot(xx-dx/2,1.2*dlineoutx,color='blue',label = 'DSMC ppc=27')
    ax.plot(xx-dx/2,1.2*dlineoutx2,color='green',label='DSMC ppc=125')
    ax.plot(xx-dx/2,1.2*dlineoutx3,color='red',label='DSMC ppc=1000')
    ax.set_xlabel('x')
    ax.set_ylabel('Density')
    ax.set_title('Density(x) @y=0, @z=0')
    plot.legend()
    plot.show()


if __name__ == '__main__':
    main()
