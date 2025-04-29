from ocelot import *
from ocelot.cpbd.elements import *
from ocelot.gui.accelerator import *
from ocelot.cpbd.track import track_debug
import matplotlib.pyplot as plt
from ocelot.cpbd.physics_proc import PhysProc

from linaccopy import *

import time
from copy import deepcopy
import numpy as np
from ocelot.cpbd.physics_proc import PhysProc




# === Define Cavity ===
cavity_voltage = 5e6  # 5 MV (can adjust)
cavity = Cavity(l=10, eid="CAV1")
cavity.voltage = cavity_voltage
cavity.freq = 1.3e9  # 1.3 GHz
cavity.phase = -20.0  # <<< Off-crest phase to introduce chirp!

# === Define Chicane (4 Dipoles) ===
# Parameters
bend_length = 1.0    # meters
bend_angle = 0.05    # radians (about 2.8 degrees)

# Bends
b1 = SBend(l=bend_length, angle=+bend_angle, eid="B1")
b2 = SBend(l=bend_length, angle=-bend_angle, eid="B2")
b3 = SBend(l=bend_length, angle=-bend_angle, eid="B3")
b4 = SBend(l=bend_length, angle=+bend_angle, eid="B4")

# Drifts between bends
dr1 = Drift(l=1.0, eid="D1")   # between B1 and B2
dr2 = Drift(l=2.0, eid="D2")   # middle drift between B2 and B3
dr3 = Drift(l=1.0, eid="D3")   # between B3 and B4

# === Define Drift after Chicane ===
drift_length = 50.0  # meters
drift_after = Drift(l=drift_length, eid="D4")

# === Build Full Toy Lattice ===
toy_lat = MagneticLattice([
    cavity,
    b1, dr1, b2, dr2, b3, dr3, b4,   # chicane section
    drift_after                      # long drift after compression
])

# Final lattice you will track on:
lat = toy_lat



##### PA1RFGUN #####
AX =-21.12   
BX =24.36 
AY =-21.16     
BY =24.42      
EMITX =1.2023504614117648e-06
EMITY =1.2023504614117648e-06 
DP =.01
energy = 0.0091
##### PA1RFGUN #####


tw0 = Twiss()
tw0.alpha_x = AX
tw0.beta_x = BX

tw0.alpha_y = AY
tw0.beta_y = BY
tw0.E = energy
tw0.emit_xn = EMITX 
tw0.emit_yn = EMITY
tw0.emit_x = EMITX   / energy * 0.511e-3
tw0.emit_y = EMITY / energy * 0.511e-3

tws = twiss(lat,tw0)

plot_opt_func(lat, tws,legend=False, grid=False, top_plot=['Dx','Dy'])

plt.show

def gauss_from_twiss_vec(emit, beta, alpha, size=1):
    phi = 2 * np.pi * np.random.rand(size)
    u = np.random.rand(size)
    a = np.sqrt(-2 * np.log(1 - u) * emit)
    x = a * np.sqrt(beta) * np.cos(phi)
    xp = -a / np.sqrt(beta) * (np.sin(phi) + alpha * np.cos(phi))
    return x, xp

nparticles = 100000
p_array = ParticleArray(n=nparticles)

# sigma_tau = 6.7e-3  # (bunch length in)
# sigma_p = 5.e-4     # relative momentum spread (Δp/p)  

sigma_tau = 4e-3  # (bunch length in )
sigma_p = 0.001     # relative momentum spread (Δp/p)  


p_array.E = energy  # GeV

# Longitudinal
p_array.rparticles[4] = np.random.randn(nparticles) * sigma_tau
p_array.rparticles[5] = np.random.randn(nparticles) * sigma_p

# Transverse: vectorized sampling from Twiss parameters
x, xp = gauss_from_twiss_vec(tw0.emit_x, tw0.beta_x, tw0.alpha_x, size=nparticles)
y, yp = gauss_from_twiss_vec(tw0.emit_y, tw0.beta_y, tw0.alpha_y, size=nparticles)

p_array.rparticles[0] = x
p_array.rparticles[1] = xp
p_array.rparticles[2] = y
p_array.rparticles[3] = yp

# Charge distribution
charge = 4.0e-9  # nC
p_array.q_array = np.ones(nparticles) * (charge / nparticles)

from ocelot.cpbd.physics_proc import PhysProc

class CustomLongitudinalWake(PhysProc):
    def __init__(self, beam_radius=0.01, bunch_length=4e-3):
        super().__init__()
        self.Z0 = 377  # Ohm
        self.a = beam_radius
        self.L = bunch_length
        self.N = None

    def apply(self, p_array, step_size):
        z = p_array.rparticles[4]
        q = p_array.q_array
        
        if self.N is None:
            self.N = len(z)
        
        z_mean = np.mean(z)
        s = z - z_mean  # centered z
        
        # Resistive Wall Wake approximation: ~ 1/sqrt(-z)
        wake = np.zeros_like(s)
        mask = s < 0  # only trailing particles feel the wake
        wake[mask] = 1.0 / np.sqrt(-s[mask] + 1e-12)  # add epsilon to avoid division by zero
        
        # Normalize and scale
        W0 = 1e3  # Scaling factor [V/C/m^1/2], can adjust
        voltage = - (self.N * np.abs(q[0]) * self.Z0 * W0 / (self.a**2)) * wake


        
        # Apply kick
        delta_E = voltage * step_size
        E = p_array.E  # GeV
        delta_dp_p = delta_E / (E * 1e9)
        
        p_array.rparticles[5] += delta_dp_p
        
        import time
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
from ocelot.cpbd.physics_proc import PhysProc

# Define the HookProcess to monitor sigma_tau
class HookProcess(PhysProc):
    def __init__(self, sigma_tau_list):
        super().__init__()
        self.sigma_tau_list = sigma_tau_list

    def apply(self, p_array, step_size):
        sigma_tau = np.std(p_array.rparticles[4])
        self.sigma_tau_list.append(sigma_tau)
        print(f"s = {self.z0:.3f} m, σ_τ = {sigma_tau:.3e} s")

# === TRACK WITHOUT WAKEFIELD ===
p_array_no_wake = deepcopy(p_array)

navi_no_wake = Navigator(lat)
navi_no_wake.unit_step = 0.5

sigma_tau_list_no_wake = []
hook_no_wake = HookProcess(sigma_tau_list_no_wake)
navi_no_wake.add_physics_proc(hook_no_wake, lat.sequence[0], lat.sequence[-1])

tws_no_wake, p_array_no_wake = track(lat, p_array_no_wake, navi_no_wake)

# === TRACK WITH CUSTOM WAKEFIELD ===
p_array_wake = deepcopy(p_array)

navi_wake = Navigator(lat)
navi_wake.unit_step = 0.5

sigma_tau_list_wake = []
hook_wake = HookProcess(sigma_tau_list_wake)
navi_wake.add_physics_proc(hook_wake, lat.sequence[0], lat.sequence[-1])

# Add custom wake
wake_proc = CustomLongitudinalWake(beam_radius=0.01, bunch_length=sigma_tau)
navi_wake.add_physics_proc(wake_proc, lat.sequence[0], lat.sequence[-1])

tws_wake, p_array_wake = track(lat, p_array_wake, navi_wake)

# === PLOT THE COMPARISON ===
s_axis_no_wake = np.linspace(0, lat.totalLen, len(sigma_tau_list_no_wake))
s_axis_wake = np.linspace(0, lat.totalLen, len(sigma_tau_list_wake))

plt.plot(s_axis_no_wake, sigma_tau_list_no_wake, label='No Wake')
plt.plot(s_axis_wake, sigma_tau_list_wake, '--', label='With Custom Wake')
plt.xlabel("s [m]")
plt.ylabel("σ_τ [s]")
plt.title("Evolution of Bunch Length σ_τ With and Without Wakefield")
plt.grid(True)
plt.legend()
plt.show()

# === PHASE SPACE PLOTS ===
def plot_phase_space(z, dp_p, title):
    plt.figure(figsize=(8,6))
    hb = plt.hexbin(z*1e4, dp_p, gridsize=100, cmap='plasma', bins='log')
    plt.colorbar(hb, label='log(Particle density)')
    plt.xlabel('z [mm]')
    plt.ylabel('Δp/p')
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.show()


# No wake
plot_phase_space(p_array_no_wake.rparticles[4], p_array_no_wake.rparticles[5], 'No Wake')

# With wake

plot_phase_space(p_array_wake.rparticles[4], p_array_wake.rparticles[5], 'With Custom Wake')


import numpy as np
import matplotlib.pyplot as plt

# Extract z positions (longitudinal coordinates) in meters
z_meters = p_array_wake.rparticles[4]  # assuming this is AFTER tracking

# Histogram to get longitudinal density profile
nbins = 2048
zmin, zmax = z_meters.min(), z_meters.max()
hist, bin_edges = np.histogram(z_meters, bins=nbins, range=(zmin, zmax))
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

# Perform FFT
fft_result = np.fft.fft(hist - np.mean(hist))  # remove DC offset
fft_freq = np.fft.fftfreq(nbins, d=(bin_centers[1] - bin_centers[0]))  # in 1/m

# Compute bunching factor magnitude
bunching = np.abs(fft_result) / np.sum(hist)

# Plot
plt.figure(figsize=(10, 6))
plt.plot(fft_freq[1:nbins//2], bunching[1:nbins//2])
plt.xlabel("Wavenumber $k$ [1/m]")
plt.ylabel("Bunching factor $b(k)$")
plt.title("Fourier-based Microbunching Analysis")
plt.grid(True)
plt.tight_layout()
plt.show()


