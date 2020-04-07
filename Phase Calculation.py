from Simulation import *
from Constants import *
import numpy as np
from copy import copy

#recoil units:
#[L] = 2/k_L
#[T] = omega_r^-1
#[f] = omega_r
#[k] = k_L/2


# class Runge_Kutta:
#
#     def __init__(self, dy):


#gamma = 1e-15

def x_th(x0, v0, t):
    x = np.zeros(3)
    x[0] = x0[0]+v0[0]*t
    x[1] = x0[1]+v0[1]*t
    x[2] = g/gamma+(x0[2]-g/gamma)*np.cosh(t*np.sqrt(gamma))+v0[2]/np.sqrt(gamma)*np.sinh(t*np.sqrt(gamma))
    return x

def v_th(x0, v0, t):
    v = np.zeros(3)
    v[0] = v0[0]
    v[1] = v0[1]
    v[2] = np.sqrt(gamma)*(x0[2]-g/gamma)*np.sinh(t*np.sqrt(gamma))+v0[2]*np.cosh(t*np.sqrt(gamma))
    return v

def x_0(x0, v0, t):
    x = np.zeros(3)
    x[0] = x0[0]+v0[0]*t
    x[1] = x0[1]+v0[1]*t
    x[2] = x0[2]+v0[2]*t-0.5*g*t**2
    return x

def v_0(x0, v0, t):
    v = np.zeros(3)
    v[0] = v0[0]
    v[1] = v0[1]
    v[2] = v0[2]-g*t
    return v

def L_th(x, v):
    return np.sum(v**2)-2*g*x[2]+gamma*x[2]**2

x_fn = x_th
v_fn = v_th
L_fn = lambda x0, v0, t: L_th(x_fn(x0, v0, t), v_fn(x0, v0, t))

class AI:

    def __init__(self, state):
        self.head = 0
        self.tail = 0
        self.states = [state]

    def add_state(self, state):
        self.tail = self.tail+1
        self.states.append(state)
        print('index = ', self.tail, 't = ', state.t, ' v = ', state.v[2],' z = ', state.x[2], ' phase = ', state.phase)

    def pop_state(self):
        self.head = self.head+1
        return self.states[self.head-1]


    def implement(self, pulse_sequence):
        last_t = 0
        for pulse in pulse_sequence:
            if pulse.t >= last_t:
                while last_t+dt < pulse.t:
                    self.free_evolution(last_t+dt)
                    last_t = last_t+dt
                self.free_evolution(pulse.t)
                self.beam_splitting(pulse)
                last_t = pulse.t
        while last_t+dt < t_det:
            self.free_evolution(last_t+dt)
            last_t = last_t+dt
        if t_det > last_t:
            self.free_evolution(t_det)

    def free_evolution(self, t):
        tail = self.tail
        while self.head <= tail:
            current_state = self.pop_state()
            self.add_state(current_state.free_evolution(t-current_state.t))

    def beam_splitting(self, pulse):
        tail = self.tail
        while self.head <= tail:
            current_state = self.pop_state()
            new_states = current_state.beam_splitting(pulse)
            for new_state in new_states:
                self.add_state(new_state)


class state:

    def __init__(self, **kwargs):
        self.x = kwargs['x'] if 'x' in kwargs else np.array([0,0,0])
        self.t = kwargs['t'] if 't' in kwargs else 0
        self.v = kwargs['v'] if 'v' in kwargs else np.array([0,0,0])
        self.amplitude = kwargs['amplitude'] if 'amplitude' in kwargs else 0
        self.phase = kwargs['phase'] if 'phase' in kwargs else 0
        self.parent = None

    def free_evolution(self, T):
        new_state = state(x = self.x, t = self.t, v = self.v, amplitude = self.amplitude, phase = self.phase)
        new_state.parent = self
        x0 = self.x
        v0 = self.v
        new_state.phase = self.phase+integrate.quad(lambda t: L_fn(x0, v0, t), 0, T)[0]
        new_state.x = x_fn(x0, v0, T)
        new_state.v = v_fn(x0, v0, T)
        new_state.t = self.t+T
        return new_state

    def beam_splitting(self, pulse):
        #Doppler shift
        delay = (2*(mirrorz-fiberz)-self.x[2])/c
        #delay = 0
        t = pulse.t+delay
        omega1 = pulse.omega1+pulse.ramp1*t-pulse.k1*self.v[2]
        omega2 = pulse.omega2+pulse.ramp2*t+pulse.k2*self.v[2]
        #If not resonant, no transition

        if abs(omega1-omega2-n*pulse.k1-n*pulse.k2)<0.1:
            undeflected_state = self.free_evolution(delay)
            undeflected_state.parent = self
            deflected_state = state(x = copy(undeflected_state.x), t = undeflected_state.t, v = copy(undeflected_state.v), amplitude = undeflected_state.amplitude, phase = undeflected_state.phase)
            deflected_state.parent = self
            deflected_state.v[2] = deflected_state.v[2]+(n*pulse.k1+n*pulse.k2)/2
            deflected_state.phase = deflected_state.phase+n*((pulse.k1+pulse.k2)*deflected_state.x[2]-(pulse.omega1-pulse.omega2)*t-pulse.phase10-0.5*pulse.ramp1*t**2+pulse.phase20+0.5*pulse.ramp2*t**2)
            print('currentz = ',deflected_state.x[2])
            print('k phase+=',n*(pulse.k1+pulse.k2)*deflected_state.x[2])
            print('ramp phase+=',n*(0.5*pulse.ramp1*t**2-0.5*pulse.ramp2*t**2))
            #deflected_state.phase = deflected_state.phase+0*n*((pulse.k1+pulse.k2)*deflected_state.x[2]+(pulse.omega1-pulse.omega2)*t+pulse.phase10-pulse.phase20)
            return undeflected_state, deflected_state
        if abs(omega2-omega1-n*pulse.k1-n*pulse.k2)<0.1:
            undeflected_state = self.free_evolution(delay)
            undeflected_state.parent = self
            deflected_state = state(x = copy(undeflected_state.x), t = undeflected_state.t, v = copy(undeflected_state.v), amplitude = undeflected_state.amplitude, phase = undeflected_state.phase)
            deflected_state.parent = self
            deflected_state.v[2] = deflected_state.v[2]-(n*pulse.k1+n*pulse.k2)/2
            deflected_state.phase = deflected_state.phase-n*((pulse.k1+pulse.k2)*deflected_state.x[2]-(pulse.omega1-pulse.omega2)*t-pulse.phase10-0.5*pulse.ramp1*t**2+pulse.phase20+0.5*pulse.ramp2*t**2)
            print('currentz = ',deflected_state.x[2])
            print('k phase-=',n*(pulse.k1+pulse.k2)*deflected_state.x[2])
            print('laser phase-=',n*(0.5*pulse.ramp1*t**2-0.5*pulse.ramp2*t**2))
            #deflected_state.phase = deflected_state.phase-0*n*((pulse.k1+pulse.k2)*deflected_state.x[2]+(pulse.omega1-pulse.omega2)*t+pulse.phase10-pulse.phase20)
            return undeflected_state, deflected_state
        return ()

    def trace_history(self):
        if self.parent == None:
            x_list = [self.x]
            t_list = [self.t]
            v_list = [self.v]
            amp_list = [self.amplitude]
            phase_list = [self.phase]
            return x_list, t_list, v_list, amp_list, phase_list
        x_list, t_list, v_list, amp_list, phase_list = self.parent.trace_history()
        x_list.append(self.x)
        t_list.append(self.t)
        v_list.append(self.v)
        amp_list.append(self.amplitude)
        phase_list.append(self.phase)
        return x_list, t_list, v_list, amp_list, phase_list

class pulse:

    def __init__(self, **kwargs):
        self.omega0 = kwargs['omega0'] if 'omega0' in kwargs else omega0
        self.omega1 = kwargs['omega1'] if 'omega1' in kwargs else 0
        self.omega2 = kwargs['omega2'] if 'omega2' in kwargs else 0
        self.ramp1 = kwargs['ramp1'] if 'ramp1' in kwargs else 0
        self.ramp2 = kwargs['ramp2'] if 'ramp2' in kwargs else 0
        self.t = kwargs['t'] if 't' in kwargs else 0
        self.phase10 = kwargs['phase10'] if 'phase10' in kwargs else 0
        self.phase20 = kwargs['phase20'] if 'phase20' in kwargs else 0
        self.k1 = (self.omega0+self.omega1)/c
        self.k2 = (self.omega0+self.omega2)/c




#AI specific constants:
n = 5
N = 125
fiberz = conversion(-1,[1,0])
mirrorz = conversion(2,[1,0])
dt = 100
g = 0
#g = conversion(g0,[1,-2])
ramp = 2*g
#ramp = 0
v0 = -5
gamma = conversion(gamma0,[0,-2])
#v0 = conversion(2.8,[1,-1])
#v0 = -5
T = conversion(0.1,[0,1])
#Initialization
init_state = state(x = np.array([0.0,0.0,0.0]), t = 0.0, v = np.array([0.0,0.0,v0]), amplitude = 1.0, phase = 0.0)
AI_MZ = AI(init_state)
pulse_sequence = []
#First Pulse
omega1 = 2*(v0+n)
omega2 = -2*(v0+n)
ramp1 = -ramp
ramp2 = ramp
t = 0
pulse1 = pulse(omega1 = omega1, omega2 = omega2, ramp1 = ramp1, ramp2 = ramp2, t = t)
pulse_sequence.append(pulse1)
#Second Pulse
omega1 = 2*(v0+n)
omega2 = -2*(v0+n)
ramp1 = -ramp
ramp2 = ramp
t = T
pulse2 = pulse(omega1 = omega1, omega2 = omega2, ramp1 = ramp1, ramp2 = ramp2, t = t)
pulse_sequence.append(pulse2)
#Bloch Oscillations
#Third Pulse
omega1 = 2*(v0+n)
omega2 = -2*(v0+n)
ramp1 = -ramp
ramp2 = ramp
t = 2*T
pulse3 = pulse(omega1 = omega1, omega2 = omega2, ramp1 = ramp1, ramp2 = ramp2, t = t)
pulse_sequence.append(pulse3)
#Fourth Pulse
0.0002598*25
-4*n**2*(T)
(-5)
-0.0194903866215-0.00649679
0.48607-0.51205
t_det = 2.01*T
4*n*ramp*T**2
(v0*T-0.5*g*T**2)*n*4
n*2*g*T**2
3741986+3612140+63410232
63280386-70504668
AI_MZ.implement(pulse_sequence)

543078656-528630093

len(AI_MZ.states)
AI_MZ.head
AI_MZ.tail
AI_MZ.states[0].x, AI_MZ.states[0].t, AI_MZ.states[0].v
AI_MZ.states[1].x, AI_MZ.states[1].t, AI_MZ.states[1].v
AI_MZ.states[2].x, AI_MZ.states[2].t, AI_MZ.states[2].v
AI_MZ.states[4].v
T
T*2
x1_list, t1_list, v1_list, amp1_list, phase1_list = AI_MZ.states[92].trace_history()
x1_list = np.array(x1_list)
t1_list = np.array(t1_list)
x2_list, t2_list, v2_list, amp2_list, phase2_list = AI_MZ.states[89].trace_history()
x2_list = np.array(x2_list)
t2_list = np.array(t2_list)
ax = plt.subplot()
ax.plot(t1_list, x1_list[:,2], 'o', t2_list, x2_list[:,2], '*')
plt.show()
AI_MZ.states[92].phase
AI_MZ.states[89].phase
AI_MZ.states[92].x[2]
AI_MZ.states[89].x[2]
(AI_MZ.states[148].x[2]-AI_MZ.states[145].x[2])*AI_MZ.states[148].v[2]
g
gamma
2*2*g*T**2*n
AI_MZ.states[148].x[2]
AI_MZ.states[144].x[2]
AI_MZ.states[148].x[2]-AI_MZ.states[145].x[2]
AI_MZ.states[149].x[2]-AI_MZ.states[144].x[2]
AI_MZ.states[144].t
AI_MZ.states[145].t
AI_MZ.states[148].t
AI_MZ.states[149].t
t_det
AI_MZ.states[100].x[2]-AI_MZ.states[97].x[2]
AI_MZ.states[101].x[2]-AI_MZ.states[96].x[2]
AI_MZ.states[100].v[2]-AI_MZ.states[97].v[2]
AI_MZ.states[100].v[2]-AI_MZ.states[97].v[2]
AI_MZ.states[100].phase-AI_MZ.states[97].phase
AI_MZ.states[101].phase-AI_MZ.states[96].phase
2*2*T**2*n*gamma*(7/12*g*T**2)
v0
v0
all_t = []
all_x = []
for current_state in AI_MZ.states:
    all_t.append(current_state.t)
    all_x.append(current_state.x)
all_t = np.array(all_t)
all_x = np.array(all_x)
ax = plt.subplot()
ax.plot(all_t, all_x[:,2], '.')
plt.show()
pulse_sequence[0].omega1
pulse_sequence[0].omega2
pulse1.k2
5 in None
