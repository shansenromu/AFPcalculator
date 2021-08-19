import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
#from scipy.interpolate import griddata
import math

def ReadB(infile,nx,ny,nz):
    arr = np.genfromtxt(infile)
    #extract x, y, z1
    xd = arr[:,0]  #[cm]
    yd = arr[:,1]
    zd = arr[:,2]
    Bx = arr[:,3] #[mT  or uT?]
    By = arr[:,4]
    Bz = arr[:,5]
    #number of steps
    #nx = 13
    #ny = 2
    #nz = 2

    #print "Point is: ", x, y, z

    xField = np.zeros((nx+1,ny+1,nz+1)) #[mT  or uT?]
    yField = np.zeros((nx+1,ny+1,nz+1))
    zField = np.zeros((nx+1,ny+1,nz+1))

    for iz in range(0,nz):
        for iy in range(0,ny):
            for ix in range(0,nx):
                #print ix, iy, iz
                #print int(ny*nx*iz+nx*iy+ix)
                xField[ix][iy][iz] = Bx[int(ny*nx*iz+nx*iy+ix)]
                yField[ix][iy][iz] = By[int(ny*nx*iz+nx*iy+ix)]
                zField[ix][iy][iz] = Bz[int(ny*nx*iz+nx*iy+ix)]
    return xd, yd, zd, xField, yField, zField

def ReadB1(infile, ny,nz,nx):    
    #read data points from text file, could be also a csv reader
    #arr = np.genfromtxt("CourseFieldMap-SCM18-070.txt")
    arr = np.genfromtxt(infile)
    #extract x, y, z1
    yd1 = arr[:,0]  #[cm]
    zd1 = arr[:,1]
    xd1 = arr[:,2]
    Bz1 = arr[:,3] #[mT  or uT?] #already aranged this from the rootfile to tsv file so this was orgianally Bz , Bx, By order
    Bx1 = arr[:,4]
    By1 = arr[:,5]


    #number of steps
    #nx = 13
    #ny = 2
    #nz = 2

    #print "Point is: ", x, y, z

    xField1 = np.zeros((nx+1,ny+1,nz+1)) #[mT  or uT?]
    yField1 = np.zeros((nx+1,ny+1,nz+1))
    zField1 = np.zeros((nx+1,ny+1,nz+1))

    for iz in range(0,nz):
        for iy in range(0,ny):
            for ix in range(0,nx):
                #print ix, iy, iz
                #print int(ny*nx*iz+nx*iy+ix)
                xField1[ix][iy][iz] = Bx1[int(ny*nx*iz+nx*iy+ix)]
                yField1[ix][iy][iz] = By1[int(ny*nx*iz+nx*iy+ix)]
                zField1[ix][iy][iz] = Bz1[int(ny*nx*iz+nx*iy+ix)]
    return xd1,yd1,zd1, xField1, yField1, zField1



def bField(x0d,y0d,z0d,xField00,yField00,zField00, x,y,z, nx,ny,nz):    
  
    xd = x0d  #[cm]
    yd = y0d
    zd = z0d
    xField = xField00
    yField = yField00
    zField = zField00

    minx = np.min(xd)
    maxx = np.max(xd)
    miny = np.min(yd)
    maxy = np.max(yd)
    minz = np.min(zd)
    maxz = np.max(zd)
    
    #distances field spans
    dx = maxx-minx
    dy = maxy-miny
    dz = maxz-minz

    #get B0 from point
    #x = 10.0
    #y = 0.0
    #z = 0.0
    point = [x,y,z]
    Bfield = [0.0,0.0,0.0]


    # Check that the point is within the defined region  
    if x>=minx and x<=maxx and \
       y>=miny and y<=maxy and \
       z>=minz and z<=maxz : 
        #Position of given point within region, normalized to the range [0,1]     
        xfraction = (x - minx) / dx     
        yfraction = (y - miny) / dy     
        zfraction = (z - minz) / dz     
        #if (invertX) { xfraction = 1 - xfraction;}     
        #if (invertY) { yfraction = 1 - yfraction;}     
        #if (invertZ) { zfraction = 1 - zfraction;}     
         
                   
        # Position of the point within the cuboid defined by the     
        # nearest surrounding tabulated points 
        modx = math.modf(xfraction*(nx-1))
        mody = math.modf(yfraction*(ny-1))
        modz = math.modf(zfraction*(nz-1))
        xlocal =  modx[0]     
        ylocal =  mody[0]    
        zlocal =  modz[0] 
        xdindex = modx[1]
        ydindex = mody[1]
        zdindex = modz[1] 

        #print xfraction, yfraction, zfraction
        #print xlocal, ylocal, zlocal
        #print (xdindex, ydindex, zdindex)
       
        #The indices of the nearest tabulated point whose coordinates are all less than those of the given point     
        xindex = int(xdindex)     
        yindex = int(ydindex)     
        zindex = int(zdindex)     
        #print (xindex, yindex, zindex)
  


        #Full 3-dimensional version    
        Bfield[0] =   \
            xField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) + \
            xField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  + \
            xField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) + \
            xField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  + \
            xField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) + \
            xField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  + \
            xField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) + \
            xField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal      
        Bfield[1] =   \
            yField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) + \
            yField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  + \
            yField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) + \
            yField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  + \
            yField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) + \
            yField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  + \
            yField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) + \
            yField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal      
        Bfield[2] =    \
            zField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) + \
            zField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  + \
            zField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) + \
            zField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  + \
            zField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) + \
            zField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  + \
            zField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) + \
            zField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal   
        
    else :
        Bfield[0] = 0.0     
        Bfield[1] = 0.0     
        Bfield[2] = 0.0
        
    #print "point, Bfield: ", point, Bfield 
    return Bfield






#Btest = bField("CourseFieldMap-SCM18-070.txt",10.0,-5.0,0.0, 13,2,2)   
#print Btest[0], Btest[1],Btest[2]


def b1Field(x1d,y1d,z1d,xField11,yField11,zField11, x,y,z, ny,nz,nx):   

    #Read in arrays of x, y, z to determine spacing
    xd1 = x1d #cm
    yd1 = y1d
    zd1 = z1d
    #Read in Feild maps xField[i][j][k]
    xField1 = xField11
    yField1 = yField11
    zField1 = zField11

  
    minx1 = np.min(xd1)
    maxx1 = np.max(xd1)
    miny1 = np.min(yd1)
    maxy1 = np.max(yd1)
    minz1 = np.min(zd1)
    maxz1 = np.max(zd1)
    
    #distances field spans
    dx1 = maxx1-minx1
    dy1 = maxy1-miny1
    dz1 = maxz1-minz1

    #get B0 from point
    #x = 10.0
    #y = 0.0
    #z = 0.0
    point = [x,y,z]
    B1field = [0.0,0.0,0.0]


    # Check that the point is within the defined region  
    if x>=minx1 and x<=maxx1 and \
       y>=miny1 and y<=maxy1 and \
       z>=minz1 and z<=maxz1 : 
        #Position of given point within region, normalized to the range [0,1]     
        xfraction = (x - minx1) / dx1     
        yfraction = (y - miny1) / dy1     
        zfraction = (z - minz1) / dz1     
        #if (invertX) { xfraction = 1 - xfraction;}     
        #if (invertY) { yfraction = 1 - yfraction;}     
        #if (invertZ) { zfraction = 1 - zfraction;}     
         
                   
        # Position of the point within the cuboid defined by the     
        # nearest surrounding tabulated points 
        modx = math.modf(xfraction*(nx-1))
        mody = math.modf(yfraction*(ny-1))
        modz = math.modf(zfraction*(nz-1))
        xlocal =  modx[0]     
        ylocal =  mody[0]    
        zlocal =  modz[0] 
        xdindex = modx[1]
        ydindex = mody[1]
        zdindex = modz[1] 

        #print xfraction, yfraction, zfraction
        #print xlocal, ylocal, zlocal
        #print (xdindex, ydindex, zdindex)
       
        #The indices of the nearest tabulated point whose coordinates are all less than those of the given point     
        xindex = int(xdindex)     
        yindex = int(ydindex)     
        zindex = int(zdindex)     
        #print (xindex, yindex, zindex)
  


        #Full 3-dimensional version    
        B1field[0] =   \
            xField1[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) + \
            xField1[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  + \
            xField1[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) + \
            xField1[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  + \
            xField1[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) + \
            xField1[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  + \
            xField1[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) + \
            xField1[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal      
        B1field[1] =   \
            yField1[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) + \
            yField1[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  + \
            yField1[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) + \
            yField1[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  + \
            yField1[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) + \
            yField1[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  + \
            yField1[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) + \
            yField1[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal      
        B1field[2] =    \
            zField1[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) + \
            zField1[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  + \
            zField1[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) + \
            zField1[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  + \
            zField1[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) + \
            zField1[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  + \
            zField1[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) + \
            zField1[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal   
        
    else :
        B1field[0] = 0.0     
        B1field[1] = 0.0     
        B1field[2] = 0.0
        
    #print "point, Bfield: ", point, Bfield 
    return B1field







#bField("CourseFieldMap-SCM18-070.txt",20.0,-5.0,0.0, 13,2,2)  
#bField("CourseFieldMap-SCM18-070.txt",10.0,4.999,0.0, 13,2,2)   
#bField("CourseFieldMap-SCM18-070.txt",20.0,4.999,0.0, 13,2,2)

