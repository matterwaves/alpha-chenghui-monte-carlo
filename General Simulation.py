from Simulation import *

#Curvature phase simulation 110317
eta = np.linspace(0.01,0.5,100)
contrast_up = (4*eta**2*(1-eta)**2)/(eta**3*(1-eta)+eta*(1-eta)**3+2*eta**2*(1-eta)**2)
contrast_down = (4*eta**2*(1-eta)**2)/(eta**4+(1-eta)**4+eta**3*(1-eta)+eta*(1-eta)**3)
ax = plt.subplot()
ax.plot(eta, contrast_up, 'b-', eta, contrast_down, '-r')
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
def Bragg_efficiency(r, w0):
    Omega_2ph_rel = np.exp(-2*r**2/w0**2)
    eta0 = eff0(Omega_2ph_rel)
    eta1 = eff1(Omega_2ph_rel)
    return eta0, eta1
#Scenario 1: no thermal expansion
def contrast_up(r, w0):
    eta0, eta1 = Bragg_efficiency(r, w0)
    return (4*eta1**2*eta0**2)/(eta1**3*eta0+eta1*eta0**3+2*eta1**2*eta0**2)
def contrast_down(r, w0):
    eta0, eta1 = Bragg_efficiency(r, w0)
    return (4*eta1**2*eta0**2)/(eta1**4+eta0**4+eta1**3*eta0+eta1*eta0**3)


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
