from Simulation import *



#Beam structure and phase
#from data to construct 2D
#Input
k = 2*np.pi/(852*10**(-9))
beam1d_x_in = np.loadtxt('beam1d_x_in.txt')
beam1d_y_in = np.loadtxt('beam1d_y_in.txt')
x_in = beam1d_x_in[450:7050,0]
signalx_in = beam1d_x_in[450:7050,1]/100
y_in = beam1d_y_in[450:7050,0]
signaly_in = beam1d_y_in[450:7050,1]/100
Ax_in = np.zeros(len(x_in))
Ay_in = np.zeros(len(y_in))
for i in range(len(x_in)):
    if signalx_in[i]>0:
        Ax_in[i] = np.sqrt(signalx_in[i])
        Ay_in[i] = np.sqrt(signaly_in[i])

phasex_in = np.zeros(len(x_in)-2)
for i in range(len(x_in)):
    if i > 0 and i < len(x_in)-1:
        phasex_in[i-1] = (Ax_in[i-1]+Ax_in[i+1]-2*Ax_in[i])/(x_in[1]-x_in[0])**2/Ax_in[i]/2/k**2
phasey_in = np.zeros(len(y_in)-2)
for i in range(len(y_in)):
    if i > 0 and i < len(y_in)-1:
        phasey_in[i-1] = (Ay_in[i-1]+Ay_in[i+1]-2*Ay_in[i])/(y_in[1]-y_in[0])**2/Ay_in[i]/2/k**2
x_in = x_in[1:-1]
Ax_in = Ax_in[1:-1]
y_in = y_in[1:-1]
Ay_in = Ay_in[1:-1]
plt.plot(x_in, Ax_in)
plt.show()
plt.plot(x_in, phasex_in)
plt.show()

#Output
beam1d_x_out = np.loadtxt('beam1d_x_out.txt')
beam1d_y_out = np.loadtxt('beam1d_y_out.txt')
x_out = beam1d_x_out[450:7050,0]
signalx_out = beam1d_x_out[450:7050,1]/100
y_out = beam1d_y_out[450:7050,0]
signaly_out = beam1d_y_out[450:7050,1]/100
Ax_out = np.zeros(len(x_out))
Ay_out = np.zeros(len(y_out))
for i in range(len(x_out)):
    if signalx_out[i]>0:
        Ax_out[i] = np.sqrt(signalx_out[i])
        Ay_out[i] = np.sqrt(signaly_out[i])

phasex_out = np.zeros(len(x_out)-2)
for i in range(len(x_out)):
    if i > 0 and i < len(x_out)-1:
        phasex_out[i-1] = (Ax_out[i-1]+Ax_out[i+1]-2*Ax_out[i])/(x_out[1]-x_out[0])**2/Ax_out[i]/2/k**2
phasey_out = np.zeros(len(y_out)-2)
for i in range(len(y_out)):
    if i > 0 and i < len(y_out)-1:
        phasey_out[i-1] = (Ay_out[i-1]+Ay_out[i+1]-2*Ay_out[i])/(y_out[1]-y_out[0])**2/Ay_out[i]/2/k**2

x_out = x_out[1:-1]
Ax_out = Ax_out[1:-1]
y_out = y_out[1:-1]
Ay_out = Ay_out[1:-1]
plt.plot(y_out, Ay_out)
plt.show()
plt.plot(y_out, phasey_out)
plt.show()
#Fit Input
def Hermite_Gaussian(x, x0, w, c0, c1, c2):
    return hermite.hermval((x-x0)*np.sqrt(2)/w,[c0,c1,c2])*np.exp(-(x-x0)**2/w**2)
p = 4
def DualGaussian(x, a1, x01, w1, a2, x02, w2):
    return a1*np.exp(-(x-x01)**2/w1**2)+a2*np.exp(-(x-x02)**2/w2**2)
def SuperGaussian(x, a, x0, w, p):
    return a*np.exp(-(x-x0)**p/w**p)

popt, pcov = curve_fit(DualGaussian, x_in, Ax_in)
popt
#popt = np.array([1,150,3600,6])
plt.plot(x_in, Ax_in, 'bo', x_in, DualGaussian(x_in, *popt), 'r-')
plt.show()
def Hermite_Gaussian2(x, x0, w, c0):
    return hermite.hermval((x-x0)*np.sqrt(2)/w,[c0,c1,c2])*np.exp(-(x-x0)**2/w**2)
c1 = -0.00180625
c2 = -0.11087475
popt, pcov = curve_fit(DualGaussian, x_in, Ax_in)
popt
plt.plot(x_in, Ax_in, 'bo', x_in, DualGaussian(x_in, *popt), 'r-')
plt.show()

popt, pcov = curve_fit(DualGaussian, -x_out, Ax_out)
popt
plt.plot(-x_out, Ax_out, 'bo', x_out, DualGaussian(x_out, *popt), 'r-')
plt.show()

c1 = 0
c2 = 0
def Hermite_Gaussian2(x, x0, w, c0):
    return hermite.hermval((x-x0)*np.sqrt(2)/w,[c0,c1,c2])*np.exp(-(x-x0)**2/w**2)
popt, pcov = curve_fit(DualGaussian, y_in, Ay_in)
popt
plt.plot(y_in, Ay_in, 'bo', y_in, DualGaussian(y_in, *popt), 'r-')
plt.show()

popt, pcov = curve_fit(Hermite_Gaussian, y_out, Ay_out)
popt
plt.plot(y_out, Ay_out, 'bo', y_out, Hermite_Gaussian(y_out, *popt), 'r-')
plt.show()
#Coalition Fit
np.shape(x_in)
x_in[6597]
x_in[0]

from __future__ import print_function
from lmfit import Parameters, minimize, fit_report
from numpy import random, linspace, pi, exp, sin, sign

def Hermite_Gaussian2(x, a, x0, w, c0, c1, c2, c3, c4):
    return a*hermite.hermval((x-x0)*np.sqrt(2)/w,[c0,c1,c2,c3,c4])*np.exp(-(x-x0)**2/w**2)

def residual(pars, x, data=None):
    vals = pars.valuesdict()
    x0in = vals['x0in']
    x0out = vals['x0out']
    w0 = vals['w0']
    z0 = vals['z0']
    c0 = vals['c0']
    c1 = vals['c1']
    c2 = vals['c2']

    zr = k*w0**2/2
    win = w0*np.sqrt(1+(0-z0)**2/zr**2)
    wout = w0*np.sqrt(1+(5-z0)**2/zr**2)
    model = np.array([w0/win*hermite.hermval((x[0]-x0in)*np.sqrt(2)/win,[c0,c1,c2])*np.exp(-(x[0]-x0in)**2/win**2),w0/wout*hermite.hermval((x[1]-x0out)*np.sqrt(2)/wout,[c0,c1,c2])*np.exp(-(x[1]-x0out)**2/wout**2)])
    if data is None:
        return model
    return (model[0]-data[0])**2+(model[1]-data[1])**2

fit_params = Parameters()
fit_params.add('x0in', value=-400)
fit_params.add('x0out', value=-100)
fit_params.add('w0', value=5000)
fit_params.add('z0', value=0)
fit_params.add('c0', value=0.8)
fit_params.add('c1', value=0.1)
fit_params.add('c2', value=0.1)


out = minimize(residual, fit_params, args=(np.array([y_in,y_out]),), kws={'data':np.array([Ay_in, Ay_out])})

print(fit_report(out))

Hermite_Gaussian2(5000,1,0,0,2000,0,1,1,1)
def Hermite_Gaussian2(x, a, x0in, x0out, w0, z0, c0, c1, c2):
    if x < 4000:
        win = w0*np.sqrt(1+(0-z0)**2/(k*w0**2)/2)**2
        return a*w0/win*hermite.hermval((x-x0in)*np.sqrt(2)/win,[c0,c1,c2])*np.exp(-(x-x0in)**2/win**2)
    else:
        wout = w0*np.sqrt(1+(6-z0)**2/(k*w0**2)/2)**2
        return a*w0/wout*hermite.hermval((x-x0out)*np.sqrt(2)/wout,[c0,c1,c2])*np.exp(-(x-x0out)**2/wout**2)
x_inout = np.concatenate((np.array([x_in,0*np.ones(len(x_in))]), np.array([x_out,6*np.ones(len(x_in))])), axis = 1)
x_inout = np.transpose(x_inout)
Ax_inout = np.concatenate((Ax_in, Ax_out), axis = 0)
np.shape(x_inout)
np.shape(Ax_inout)

popt, pcov = curve_fit(Hermite_Gaussian2, x_inout, Ax_inout)
popt
plt.plot(x_in, Ax_in, 'bo', x_in, Hermite_Gaussian(x_in, *popt), 'r-')
plt.show()

#Monte Carlo System


#Beam profile 1d
#Input
beam1d_x_in = np.loadtxt('beam1d_x_in_Bloch.txt')
beam1d_y_in = np.loadtxt('beam1d_y_in_Bloch.txt')
x_in = beam1d_x_in[450:7050,0]/10**6
signalx_in = beam1d_x_in[450:7050,1]/100
y_in = beam1d_y_in[450:7050,0]/10**6
signaly_in = beam1d_y_in[450:7050,1]/100
plt.plot(x_in, signalx_in)
plt.show()
Ax_in = np.zeros(len(x_in))
Ay_in = np.zeros(len(y_in))
for i in range(len(x_in)):
    if signalx_in[i]>0:
        Ax_in[i] = np.sqrt(signalx_in[i])
        Ay_in[i] = np.sqrt(signaly_in[i])
popt, pcov = curve_fit(Hermite_Gaussian, x_in, Ax_in)
popt
plt.plot(x_in, Ax_in, 'bo', x_in, Hermite_Gaussian(x_in, *popt), 'r-')
plt.show()
popt, pcov = curve_fit(Hermite_Gaussian, y_in, Ay_in)
popt
plt.plot(y_in, Ay_in, 'bo', y_in, Hermite_Gaussian(y_in, *popt), 'r-')
plt.show()

phasex_in = np.zeros(len(x_in)-2)
for i in range(len(x_in)):
    if i > 0 and i < len(x_in)-1:
        phasex_in[i-1] = (Ax_in[i-1]+Ax_in[i+1]-2*Ax_in[i])/(x_in[1]-x_in[0])**2/Ax_in[i]/2/k**2
phasey_in = np.zeros(len(y_in)-2)
for i in range(len(y_in)):
    if i > 0 and i < len(y_in)-1:
        phasey_in[i-1] = (Ay_in[i-1]+Ay_in[i+1]-2*Ay_in[i])/(y_in[1]-y_in[0])**2/Ay_in[i]/2/k**2
x_in = x_in[1:-1]
Ax_in = Ax_in[1:-1]
y_in = y_in[1:-1]
Ay_in = Ay_in[1:-1]
plt.plot(x_in, Ax_in)
plt.show()
plt.plot(x_in, phasex_in)
plt.show()

#Output
beam1d_x_out = np.loadtxt('beam1d_x_out_Bloch.txt')
beam1d_y_out = np.loadtxt('beam1d_y_out_Bloch.txt')
x_out = beam1d_x_out[450:7050,0]/10**6
signalx_out = beam1d_x_out[450:7050,1]/100
y_out = beam1d_y_out[450:7050,0]/10**6
signaly_out = beam1d_y_out[450:7050,1]/100
Ax_out = np.zeros(len(x_out))
Ay_out = np.zeros(len(y_out))
for i in range(len(x_out)):
    if signalx_out[i]>0:
        Ax_out[i] = np.sqrt(signalx_out[i])
        Ay_out[i] = np.sqrt(signaly_out[i])

phasex_out = np.zeros(len(x_out)-2)
for i in range(len(x_out)):
    if i > 0 and i < len(x_out)-1:
        phasex_out[i-1] = (Ax_out[i-1]+Ax_out[i+1]-2*Ax_out[i])/(x_out[1]-x_out[0])**2/Ax_out[i]/2/k**2
phasey_out = np.zeros(len(y_out)-2)
for i in range(len(y_out)):
    if i > 0 and i < len(y_out)-1:
        phasey_out[i-1] = (Ay_out[i-1]+Ay_out[i+1]-2*Ay_out[i])/(y_out[1]-y_out[0])**2/Ay_out[i]/2/k**2

x_out = x_out[1:-1]
Ax_out = Ax_out[1:-1]
y_out = y_out[1:-1]
Ay_out = Ay_out[1:-1]

initial_guess = (-1.56782888e-04,   2.24967698e-03,   1.15792394e+00,
         4.86633574e-02,   8.78800292e-02)
popt, pcov = curve_fit(Hermite_Gaussian, x_out, Ax_out, p0 = initial_guess)
popt
plt.plot(x_out, Ax_out, 'b.', x_out, Hermite_Gaussian(x_out, *popt), 'r-')
plt.show()
initial_guess = (1.12856854e-05,   2.31617575e-03,   1.17815218e+00,
        -2.28039203e-02,   9.43527217e-02)
popt, pcov = curve_fit(Hermite_Gaussian, y_out, Ay_out, p0 = initial_guess)
popt
plt.plot(y_out, Ay_out, 'bo', y_out, Hermite_Gaussian(y_out, *popt), 'r-')
plt.show()
plt.plot(y_out, Ay_out)
plt.show()
plt.plot(y_out, phasey_out)
plt.show()

#Fit in one
#X
k
def Hermite_Gaussian_prop(xz,x0,w0,z0,c0,c1,c2):
    x = xz[0]
    z = xz[1]
    zr = k*w0**2/2
    w = w0*np.sqrt(1+(z-z0)**2/zr**2)
    return hermite.hermval((x-x0)*np.sqrt(2)/w,[c0,c1,c2])*np.exp(-(x-x0)**2/w**2)

x1 = x_in+1.56782888e-04
x2 = x_out-1.34070820e-04
xz = [np.concatenate([x1,x2]),np.concatenate([60*0.0254*np.ones(len(x_in)),148*0.0254*np.ones(len(x_out))])]
initial_guess = (0,   2.24967698e-03,   0, 1.15792394e+00,
         4.86633574e-02,   8.78800292e-02)
popt, pcov = curve_fit(Hermite_Gaussian_prop, xz, np.concatenate([Ax_in,Ax_out]), p0 = initial_guess)
popt
w0 = 0.00227517835
z0 = 1.14565545
z1 = 60*0.0254
zr = k*w0**2/2
w1 = w0*np.sqrt(1+(z1-z0)**2/zr**2)
w1
z2 = 148*0.0254
w2 = w0*np.sqrt(1+(z2-z0)**2/zr**2)
w2
plt.plot(x1, Ax_in, 'b-', x1, Hermite_Gaussian(x1, 0, w1, *popt[3:6]), 'r-')
plt.show()

plt.plot(x2, Ax_out, 'b-', x2, Hermite_Gaussian(x2, 0, w2, *popt[3:6]), 'r-')
plt.show()

#Y
y1 = y_in-1.12856854e-05
y2 = y_out-1.09024580e-04
yz = [np.concatenate([y1,y2]),np.concatenate([60*0.0254*np.ones(len(x_in)),148*0.0254*np.ones(len(x_out))])]
initial_guess = (0,   2.31617575e-03,   1, 1.17815218e+00,
        -2.28039203e-02,   9.43527217e-02)
popt, pcov = curve_fit(Hermite_Gaussian_prop, yz, np.concatenate([Ay_in,Ay_out]), p0 = initial_guess)
popt
w0 = 0.00221754139
z0 = 7.30209408
z1 = 60*0.0254
zr = k*w0**2/2
w1 = w0*np.sqrt(1+(z1-z0)**2/zr**2)
w1
z2 = 148*0.0254
w2 = w0*np.sqrt(1+(z2-z0)**2/zr**2)
w2
plt.plot(y1, Ay_in, 'b-', y1, Hermite_Gaussian(y1, 0, w1, *popt[3:6]), 'r-')
plt.show()

plt.plot(y2, Ay_out, 'b-', y2, Hermite_Gaussian(y2, 0, w2, *popt[3:6]), 'r-')
plt.show()
w1
w2
#Beam profile 2d
#Input
k = 2*np.pi/(852*10**(-9))
beam2d_in = np.genfromtxt('beam2d_in.txt',delimiter=';')
np.shape(beam2d_in)
Axy_in = np.zeros([np.shape(beam2d_in)[0]//4,np.shape(beam2d_in)[1]//4])
np.shape(Axy_in)
for i in range(np.shape(Axy_in)[0]):
    for j in range(np.shape(Axy_in)[1]):
        centerx = i*4
        centery = j*4
        Axy_in[i,j] = np.sqrt(np.mean(beam2d_in[centerx:(centerx+4),centery:(centery+4)]))

np.shape(Axy_in)
x = np.linspace(-680*(6.45*10**-6),676*(6.45*10**-6),340)
y = np.linspace(-512*(6.45*10**-6),508*(6.45*10**-6),256)
dx = x[1]-x[0]
xx, yy = np.meshgrid(x,y)
np.shape(xx)

#fig = plt.figure()
#ax = fig.gca(projection='3d')
#surf = ax.plot_surface(xx, yy, Axy_in, cmap=cm.coolwarm, linewidth=0, antialiased=False)
#plt.show()

Axy_in_fft = np.fft.fft2(Axy_in)
Axy_in_PSD = np.log(abs(Axy_in_fft)**2)
fig, ax = plt.subplots(1, 1)
ax.imshow(Axy_in_PSD, cmap=plt.cm.jet, origin='bottom')
plt.show()
Axy_in_fft[1:10,1:10]
Axy_in_PSD[1:10,1:10]


def Gaussian2D(x, a, x0, y0, w):
    return a*np.exp(-((x[0]-x0)**2+(x[1]-y0)**2)/w**2)
l = 0
xy = np.zeros([np.shape(Axy_in)[0]*np.shape(Axy_in)[1],2])
dataxy = np.zeros(np.shape(Axy_in)[0]*np.shape(Axy_in)[1])
for i in range(np.shape(Axy_in)[0]):
    for j in range(np.shape(Axy_in)[1]):
        xy[l] = [xx[i,j], yy[i,j]]
        dataxy[l] = Axy_in[i,j]
        l = l+1
xy = np.array(xy)
initial_guess = (10, 0.00039, -7.74E-5, 0.0032)
popt, pcov = curve_fit(Gaussian2D, xy.transpose(), dataxy, p0=initial_guess)
Axy_in_amp = popt[0]
Axy_in_center = popt[1:2]
Axy_in_w = popt[3]
fig, ax = plt.subplots(1, 1)
ax.imshow(Axy_in, cmap=plt.cm.jet, origin='bottom')
ax.contour(Gaussian2D([xx,yy],*popt).reshape(np.shape(Axy_in)))
plt.show()


np.shape(xx)
dx
phase2d_in = np.zeros([np.shape(Axy_in)[0]-2,np.shape(Axy_in)[1]-2])
np.shape(phase2d_in)
for i in range(np.shape(phase2d_in)[0]):
    for j in range(np.shape(phase2d_in)[1]):
        phase2d_in[i,j] = ((Axy_in[i,j+1]+Axy_in[i+2,j+1]-2*Axy_in[i+1,j+1])/dx**2+(Axy_in[i+1,j]+Axy_in[i+1,j+2]-2*Axy_in[i+1,j+1])/dx**2)/Axy_in[i+1,j+1]/2/k**2

weights = Gaussian2D([xx,yy],1,popt[1],popt[2],0.001)
#fig, ax = plt.subplots(1, 1)
#ax.imshow(weights, cmap=plt.cm.jet, origin='bottom')
#plt.show()

np.average(phase2d_in, weights = weights[1:-1,1:-1])

#Output
beam2d_out = np.loadtxt('beam2d_out.txt')
np.shape(beam2d_out)
Axy_out = np.zeros([np.shape(beam2d_out)[0]//4,np.shape(beam2d_out)[1]//4])
np.shape(Axy_out)
for i in range(np.shape(Axy_out)[0]):
    for j in range(np.shape(Axy_out)[1]):
        centerx = i*4
        centery = j*4
        Axy_out[i,j] = np.sqrt(np.mean(beam2d_out[centerx:(centerx+4),centery:(centery+4)]))

np.shape(Axy_out)
x = np.linspace(-680*(6.45*10**-6),676*(6.45*10**-6),340)
y = np.linspace(-512*(6.45*10**-6),508*(6.45*10**-6),256)
dx = x[1]-x[0]
xx, yy = np.meshgrid(x,y)
np.shape(xx)

#fig = plt.figure()
#ax = fig.gca(projection='3d')
#surf = ax.plot_surface(xx, yy, Axy_out, cmap=cm.coolwarm, linewidth=0, antialiased=False)
#plt.show()

l = 0
xy = np.zeros([np.shape(Axy_out)[0]*np.shape(Axy_out)[1],2])
dataxy = np.zeros(np.shape(Axy_out)[0]*np.shape(Axy_out)[1])
for i in range(np.shape(Axy_out)[0]):
    for j in range(np.shape(Axy_out)[1]):
        xy[l] = [xx[i,j], yy[i,j]]
        dataxy[l] = Axy_out[i,j]
        l = l+1
xy = np.array(xy)
initial_guess = (10, 0.00039, -7.74E-5, 0.0032)
popt, pcov = curve_fit(Gaussian2D, xy.transpose(), dataxy, p0=initial_guess)
Axy_out_amp = popt[0]
Axy_out_center = popt[1:2]
Axy_out_w = popt[3]

fig, ax = plt.subplots(1, 1)
ax.imshow(Axy_out, cmap=plt.cm.jet, origin='bottom')
ax.contour(Gaussian2D([xx,yy],*popt).reshape(np.shape(Axy_out)))
plt.show()


dx
phase2d_out = np.zeros([np.shape(Axy_out)[0]-2,np.shape(Axy_out)[1]-2])
np.shape(phase2d_out)
for i in range(np.shape(phase2d_out)[0]):
    for j in range(np.shape(phase2d_out)[1]):
        phase2d_out[i,j] = ((Axy_out[i,j+1]+Axy_out[i+2,j+1]-2*Axy_out[i+1,j+1])/dx**2+(Axy_out[i+1,j]+Axy_out[i+1,j+2]-2*Axy_out[i+1,j+1])/dx**2)/Axy_out[i+1,j+1]/2/k**2

weights = Gaussian2D([xx,yy],1,popt[1],popt[2],0.001)
#fig, ax = plt.subplots(1, 1)
#ax.imshow(weights, cmap=plt.cm.jet, origin='bottom')
#plt.show()

np.average(phase2d_out, weights = weights[1:-1,1:-1])

#Pure Gaussian Test
AG = Gaussian2D([xx,yy],1,0,0,0.003)
weights = Gaussian2D([xx,yy],1,0,0,0.002*np.sqrt(2))
phaseG = np.zeros([np.shape(AG)[0]-2,np.shape(AG)[1]-2])
np.shape(phaseG)
for i in range(np.shape(phaseG)[0]):
    for j in range(np.shape(phaseG)[1]):
        phaseG[i,j] = ((AG[i,j+1]+AG[i+2,j+1]-2*AG[i+1,j+1])/dx**2+(AG[i+1,j]+AG[i+1,j+2]-2*AG[i+1,j+1])/dx**2)/AG[i+1,j+1]/2/k**2
np.average(phaseG, weights = weights[1:-1,1:-1])
np.shape(phaseG)
-2/0.003**2/k**2

Axy_up = Axy_in/Axy_in_amp
#Axy_up = AG
phasexy_up = phase2d_in
#phasexy_up = phaseG
Axy_down = Axy_out/Axy_out_amp
#Axy_down = AG
phasexy_down = phase2d_out
#phasexy_down = phaseG

np.shape(phasexy_up)

def beam_up(x, y, peak):
    i = int((y+512*(6.45*10**-6))/dx)
    j = int((x+680*(6.45*10**-6))/dx)
    if i < 1 or j < 1 or i >= np.shape(phasexy_up)[0] or j >= np.shape(phasexy_up)[1]:
        return 0,0
    #print(x, y, i, j)
    amp0 = Axy_up[i,j]+(Axy_up[i,j+1]-Axy_up[i,j])/dx*(x+680*(6.45*10**-6)-j*dx)
    amp1 = Axy_up[i+1,j]+(Axy_up[i+1,j+1]-Axy_up[i+1,j])/dx*(x+680*(6.45*10**-6)-j*dx)
    amplitude = peak*(amp0+(amp1-amp0)/dx*(y+512*(6.45*10**-6)-i*dx))
    i = i-1
    j = j-1
    phi0 = phasexy_up[i,j]+(phasexy_up[i,j+1]-phasexy_up[i,j])/dx*(x+680*(6.45*10**-6)-j*dx)
    phi1 = phasexy_up[i+1,j]+(phasexy_up[i+1,j+1]-phasexy_up[i+1,j])/dx*(x+680*(6.45*10**-6)-j*dx)
    phase = phi0+(phi1-phi0)/dx*(y+512*(6.45*10**-6)-i*dx)
    return amplitude, phase

def beam_down(x, y, peak):
    i = int((y+512*(6.45*10**-6))/dx)
    j = int((x+680*(6.45*10**-6))/dx)
    if i < 1 or j < 1 or i >= np.shape(phasexy_up)[0] or j >= np.shape(phasexy_up)[1]:
        return 0,0
    amp0 = Axy_down[i,j]+(Axy_down[i,j+1]-Axy_down[i,j])/dx*(x+680*(6.45*10**-6)-j*dx)
    amp1 = Axy_down[i+1,j]+(Axy_down[i+1,j+1]-Axy_down[i+1,j])/dx*(x+680*(6.45*10**-6)-j*dx)
    amplitude = peak*(amp0+(amp1-amp0)/dx*(y+512*(6.45*10**-6)-i*dx))
    i = i-1
    j = j-1
    phi0 = phasexy_down[i,j]+(phasexy_down[i,j+1]-phasexy_down[i,j])/dx*(x+680*(6.45*10**-6)-j*dx)
    phi1 = phasexy_down[i+1,j]+(phasexy_down[i+1,j+1]-phasexy_down[i+1,j])/dx*(x+680*(6.45*10**-6)-j*dx)
    phase = phi0+(phi1-phi0)/dx*(y+512*(6.45*10**-6)-i*dx)
    return amplitude, phase

import numpy as np
beam2d_in = np.loadtxt('beam2d_in.txt')

A = np.fromfile('beam2d_in.raw', dtype='int16', sep="")
1360*1024
np.shape(A)
1392667-1392640
A_info = A[:27]
np.shape(A_info)
A_info
A = A[27:]
A=A.reshape([1024,1360])
plt.imshow(A)
plt.show()


#Read Bragg power efficiency function
Bragg_data_table = np.loadtxt('Bragg_power.txt')
intensity = Bragg_data_table[:,0]
population0 = Bragg_data_table[:,1]
population1 = Bragg_data_table[:,2]
eff0 = interpolate.interp1d(intensity, population0)
eff1 = interpolate.interp1d(intensity, population1)
x = np.linspace(0, 1.8, 500)
y = eff0(x)
plt.plot(intensity, population0, '+', x, y, '-')
plt.show()
x = np.linspace(0, 1.8, 500)
y = eff1(x)
plt.plot(intensity, population1, '+', x, y, '-')
plt.show()

def contrast_up(intensity):
    if intensity < 0.1 or intensity > 1.25:
        return 0
    eta0, eta1 = eff0(intensity), eff1(intensity)
    return (4*eta1**2*eta0**2)/(eta1**3*eta0+eta1*eta0**3+2*eta1**2*eta0**2)
def contrast_down(intensity):
    if intensity < 0.1 or intensity > 1.25:
        return 0
    eta0, eta1 = eff0(intensity), eff1(intensity)
    return (4*eta1**2*eta0**2)/(eta1**4+eta0**4+eta1**3*eta0+eta1*eta0**3)


#Parameters
sigma_x = 0.002
sigma_v = 0.0035
kB = 1.38064852e-23
M = 2.20694650e-25
Temperature = sigma_v**2*M/kB
size = 100000
xy = np.zeros([size,2])
phase = np.zeros([size,2])
weight = np.zeros([size,2])
Bragg_peak_rel = 1.1
t0 = 1.24-1.12
T = 0.080
Tprime = 0.020
xc = np.array(Axy_in_center)
vc = np.array([0, 0])
#Bloch_peak_rel =
#Monte Carlo body
np.random.seed(255349)
for i in range(size):
    x0 = np.random.normal(0, sigma_x, 2)+xc
    vT = np.random.normal(0, sigma_v, 2)+vc
    current_xy = x0
    #First pulse
    amp_up, phase_up = beam_up(current_xy[0],current_xy[1],np.sqrt(Bragg_peak_rel))
    amp_down, phase_down = beam_down(current_xy[0],current_xy[1],np.sqrt(Bragg_peak_rel))
    weight[i,0] = contrast_up(amp_up*amp_down)
    weight[i,1] = contrast_down(amp_up*amp_down)
    phase[i,0] = phase[i,0]+phase_up
    phase[i,1] = phase[i,1]+phase_down
    xy[i,:] = current_xy

np.shape(phase)
np.average(phase[:,0],weights = weight[:,0])
np.average(phase[:,1],weights = weight[:,0])
np.mean(phase[:,1])
np.average(phase[:,0],weights = weight[:,0])
np.average(phase[:,1],weights = weight[:,0])


np.average(phase[:,0])
np.average(phase[:,1])

phase[0:10,0]
xy[0:10]
x_in
-18694/10**12

ax = plt.subplot()
ax.plot(xy[:,0],xy[:,1],'bo')
plt.show()

phase[1:20,0]

def contrast_up(r, w0):
    eta0, eta1 = Bragg_efficiency(r, w0)
    return (4*eta1**2*eta0**2)/(eta1**3*eta0+eta1*eta0**3+2*eta1**2*eta0**2)
def contrast_down(r, w0):
    eta0, eta1 = Bragg_efficiency(r, w0)
    return (4*eta1**2*eta0**2)/(eta1**4+eta0**4+eta1**3*eta0+eta1*eta0**3)




r = np.sqrt(xy[:,0]**2+xy[:,1]**2)
plt.plot(xy[:,0],xy[:,1],'o')
plt.show()

plt.hist(xy[:,0])
plt.show()
mu, sigma = norm.fit(xy[:,0])
sigma



#Curvature phase simulation 110317
eta = np.linspace(0.01,0.5,100)

ax = plt.subplot()
ax.plot(eta, contrast_up, 'b-', eta, contrast_down, '-r')
plt.show()



w0 = 3.2
r = np.linspace(0, 6, 50)
ax = plt.subplot()
ax.plot(r, contrast_up(r, w0), 'b-', r, contrast_down(r, w0), '-r')
plt.show()
integrate.quad(lambda rr: rr*np.exp(-rr**2/2*2**2)*contrast_up(rr,w0), 0, 5)[0]
integrate.quad(lambda rr: np.exp(-rr**2/2*2**2)*contrast_up(rr,w0), 0, 5)
0.2090/0.5826

integrate.quad(lambda rr: rr*np.exp(-rr**2/2*2**2)*contrast_down(rr,w0), 0, 5)
integrate.quad(lambda rr: np.exp(-rr**2/2*2**2)*contrast_down(rr,w0), 0, 5)
0.1569/0.5123

#Bloch efficiency
def Bloch_efficiency(Bloch_Omega):




t = np.linspace(0.1,1,1000)
deltap = np.zeros(len(t))
r = np.zeros(len(t))
sigma = 0.002
sigma = 0.0002
w = 0.0036
sigma_eff = sigma*w/np.sqrt(w**2-2*sigma**2)
sigma
sigma_eff
ks = 1/(0.00002)
k = 2*np.pi/(852*10**(-9))
k = 2*np.pi/(780*10**(-9))
ks
for i in range(len(t)):
    rc = w*np.sqrt(-np.log(t[i]))
    r[i] = rc
    deltap[i] = -(sigma_eff/sigma)**2*(np.exp(-sigma_eff**2*ks**2/2)-np.exp(-rc**2/2/sigma_eff**2)*np.real((np.exp(-rc*ks*1j)*special.erfcx(rc/np.sqrt(2)/sigma_eff+sigma_eff*ks/np.sqrt(2)*1j)+np.exp(rc*ks*1j)*special.erfcx(rc/np.sqrt(2)/sigma_eff-sigma_eff*ks/np.sqrt(2)*1j))))*special.erf(rc/np.sqrt(2)/sigma_eff)*ks**2/4/k**2*0.1


special.erf(rc/np.sqrt(2)/sigma_eff)
(special.erf(rc/np.sqrt(2)/sigma_eff+sigma_eff*ks/np.sqrt(2)*1j)+special.erf(rc/np.sqrt(2)/sigma_eff-sigma_eff*ks/np.sqrt(2)*1j))
sigma_eff*ks
abs(special.erfcx(10+1000j))
ks
sigma_eff*ks
rc = 0.001
sigma_eff = 0.001
ks = 1/0.0001
np.exp(-sigma_eff**2*ks**2/2)
np.exp(-rc**2/2/sigma_eff**2)
np.exp(-rc**2/2/sigma_eff**2)
np.real((np.exp(rc*ks*1j)*special.erfcx(rc/np.sqrt(2)/sigma_eff+sigma_eff*ks/np.sqrt(2)*1j)+np.exp(-rc*ks*1j)*special.erfcx(rc/np.sqrt(2)/sigma_eff-sigma_eff*ks/np.sqrt(2)*1j)))


ks = 2*np.pi/0.0005
k = 2*np.pi/(852*10**(-9))
sigma = 0.002
#sigma = 0.0002
w = 0.0036
sigma_eff = sigma*w/np.sqrt(w**2-2*sigma**2)
rc = 0.001
2*np.exp(-sigma_eff**2*ks**2/2)-np.exp(-rc**2/2/sigma_eff**2)*np.real((np.exp(-rc*ks*1j)*special.erfcx(rc/np.sqrt(2)/sigma_eff+sigma_eff*ks/np.sqrt(2)*1j)+np.exp(rc*ks*1j)*special.erfcx(rc/np.sqrt(2)/sigma_eff-sigma_eff*ks/np.sqrt(2)*1j)))
-(sigma_eff/sigma)**2*(np.exp(-sigma_eff**2*ks**2/2)-np.exp(-rc**2/2/sigma_eff**2)*np.real((np.exp(-rc*ks*1j)*special.erfcx(rc/np.sqrt(2)/sigma_eff+sigma_eff*ks/np.sqrt(2)*1j)+np.exp(rc*ks*1j)*special.erfcx(rc/np.sqrt(2)/sigma_eff-sigma_eff*ks/np.sqrt(2)*1j))))*special.erf(rc/np.sqrt(2)/sigma_eff)*ks**2/4/k**2*0.1

k = 2*np.pi/(780*10**(-9))
sigma = 0.0002
w = 0.0036
sigma_eff = sigma*w/np.sqrt(w**2-2*sigma**2)
2*np.exp(-sigma_eff**2*ks**2/2)-np.exp(-rc**2/2/sigma_eff**2)*np.real((np.exp(-rc*ks*1j)*special.erfcx(rc/np.sqrt(2)/sigma_eff+sigma_eff*ks/np.sqrt(2)*1j)+np.exp(rc*ks*1j)*special.erfcx(rc/np.sqrt(2)/sigma_eff-sigma_eff*ks/np.sqrt(2)*1j)))
-(sigma_eff/sigma)**2*(np.exp(-sigma_eff**2*ks**2/2)-np.exp(-rc**2/2/sigma_eff**2)*np.real((np.exp(-rc*ks*1j)*special.erfcx(rc/np.sqrt(2)/sigma_eff+sigma_eff*ks/np.sqrt(2)*1j)+np.exp(rc*ks*1j)*special.erfcx(rc/np.sqrt(2)/sigma_eff-sigma_eff*ks/np.sqrt(2)*1j))))*special.erf(rc/np.sqrt(2)/sigma_eff)*ks**2/4/k**2*0.1




sigma_eff
rc
ks
sigma_eff*ks
1/2/sigma_eff**2
plt.plot(t, deltap)
plt.show()
deltap[980]
special.erf(1+10j)+special.erf(1-10j)
2-special.erfc(1+10j)-special.erfc(1-10j)
2*np.exp(-100/2)-np.exp(-1/2)*(np.exp(-10j)*special.erfcx((1+10j)/np.sqrt(2))+np.exp(10j)*special.erfcx((1-10j)/np.sqrt(2)))
#Scenario 2: with thermal expansion


#Delete Modules

for module in ['Simulation']:
    delete_module(module)
