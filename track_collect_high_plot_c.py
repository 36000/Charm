'''
high level variables (16):'''

high_var = ['jet_pt', 'jet_eta',
'track_2_d0_significance', 'track_3_d0_significance',
'track_2_z0_significance', 'track_3_z0_significance',
'n_tracks_over_d0_threshold', 'jet_prob', 'jet_width_eta', 'jet_width_phi',
'vertex_significance', 'n_secondary_vertices', 'n_secondary_vertex_tracks',
'delta_r_vertex', 'vertex_mass', 'vertex_energy_fraction']

'''mid level variables (28) * max tracks (15):

[D0, Z0, PHI, THETA, QOVERP,
D0D0,
Z0D0, Z0Z0,
PHID0, PHIZ0, PHIPHI,
THETAD0, THETAZ0, THETAPHI, THETATHETA,
QOVERPD0, QOVERPZ0, QOVERPPHI, QOVERPTHETA, QOVERPQOVERP,
track_weight,
mass, displacement, delta_eta_jet, delta_phi_jet,
displacement_significance, n_tracks, energy_fraction]

flavor (y):
signal --> y == 5
'''
import numpy as np
import matplotlib.pyplot as plt

histo_sig = np.genfromtxt("histo_sig_collector.csv", delimiter=',')
histo_c = np.genfromtxt("histo_c_collector.csv", delimiter=',')
histo_bg = np.genfromtxt("histo_bg_collector.csv", delimiter=',')
bino = np.genfromtxt("bin_collector.csv", delimiter=',')

for i in range(histo_sig.shape[0]):
    pltSig = plt.plot(bino[i,:][:-1], histo_sig[i,:], drawstyle='steps-post', color='blue', label='bottom')
    pltC = plt.plot(bino[i,:][:-1], histo_c[i,:], drawstyle='steps-post', color='green', label='charm')
    pltBg = plt.plot(bino[i,:][:-1], histo_bg[i,:], drawstyle='steps-post', color='red', label='bg')
    plt.title(high_var[i])
    plt.legend(loc='upper right')
    plt.savefig('h'+str(i)+'.png')
    plt.clf()
