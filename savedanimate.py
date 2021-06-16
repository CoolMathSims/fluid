from PIL import Image
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import PillowWriter
import numpy as np
import os
datu = np.loadtxt('u.dat')
datv = np.loadtxt('v.dat')
datT = np.loadtxt('eta.dat')
top = np.loadtxt('topo.dat')
top = top[1:-1, 1:-1]
x = np.linspace(0,100,100)
y = np.linspace(0,50,50)
X,Y=np.meshgrid(x,y)
levels = np.arange(0,2,0.05)

#fig,ax = plt.figure(figsize=(10,6))

images = []
#h = ax.plot_wireframe(X,Y,height,color='black')
plot = None
flow = None
topo = None
h = None
n=0
Rey = []
dim = len(datu)/100

while n<dim:
    fig,ax = plt.subplots(figsize=(10,6))
    ind = np.arange(0,50) + n*100
    Uuse = []
    Vuse = []
    Tuse = []
    for i in ind:
        u = datu[i]
        v = datv[i]
        t = datT[i]
        Uuse.append(u)
        Vuse.append(v)
        Tuse.append(t)
    Uuse = np.array(Uuse)
    Vuse = np.array(Vuse)
    Tuse = np.array(Tuse)
    incident = Uuse[:,0]
    #if flow:
    #    ax.collections=[]
    #    ax.collections.remove(flow)

    #flow = ax.contour(X,Y,Tuse,levels)
    #flow = ax.streamplot(X,Y,Uuse,Vuse)
    #topo = ax.contour(X,Y,top)
    #turb = ax.imshow(Tuse)
    Z=Tuse
    fig = plt.figure(figsize=(20,10))
    ax = plt.axes(projection='3d')
    ax.set_zlim(-1,1)
    flow = ax.plot_wireframe(X, Y, -Z, color='blue',alpha=0.5)
    im = plt.savefig('plot'+str(n)+'.png')
    images.append('plot'+str(n)+'.png')
    plt.close()


    #print(n)
    n=n+2

frames = []
for i in images:
    new_frame = Image.open(i)
    frames.append(new_frame)
    os.remove(i)

frames[0].save('animatedWire.gif', format='GIF', append_images=frames[1:],save_all=True, loop=0)
