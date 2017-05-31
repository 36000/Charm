'''
high level variables (16):

[jet_pt, jet_eta,
track_2_d0_significance, track_3_d0_significance,
track_2_z0_significance, track_3_z0_significance,
n_tracks_over_d0_threshold, jet_prob, jet_width_eta, jet_width_phi,
vertex_significance, n_secondary_vertices, n_secondary_vertex_tracks,
delta_r_vertex, vertex_mass, vertex_energy_fraction]

mid level variables (28) * max tracks (15):

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

_binning = {0:[0,200], 2:[0,100], 3:[0,50], 4:[0,100], 5:[0,100], 6:[0,10], 7:[0,0.4], 10:[0,1000], 14:[0,10]}

import numpy as np

filepath = '/phys/groups/tev/scratch4/users/kaifulam/dguest/gjj-pheno/v1/julian/raw_data'

high_sig_collect = []
high_c_collect = []
high_bg_collect = []
y_collect = []

#ONLY HIGH LEVEL VARIABLE

for i in range(10000):
    if i % 1000 == 0:
        print '@line number k = ', i

    high = np.load(filepath+'/saved_batches_test/clean_dijet_high_'+str(i)+'.npy')
    #mid = np.load(filepath+'/saved_batches_test/clean_dijet_mid_'+str(i)+'.npy')
    y = np.load(filepath+'/saved_batches_test/clean_dijet_y_'+str(i)+'.npy')
    
    # y==5 is b jet signal
    if y == 5:
        high_sig_collect.append(high)
    elif y == 4:
        high_c_collect.append(high)
    else:
        high_bg_collect.append(high)

    y_collect.append(y)

high_sig_collect = np.asarray(high_sig_collect)
high_c_collect = np.asarray(high_c_collect)
high_bg_collect = np.asarray(high_bg_collect)
y_collect = np.asarray(y_collect)

#print 'mid', mid
#print 'y_collect', y_collect

histo_sig_collector = []
histo_c_collector = []
histo_bg_collector = []
bin_collector = []

#histogram for each variable k
for k in range(high_sig_collect.shape[2]):

    var_sig = high_sig_collect[:,:,k]
    var_c = high_c_collect[:,:,k]
    var_bg = high_bg_collect[:,:,k]

    #remove nan for plotting
    sig = var_sig[~np.isnan(var_sig)]
    c = var_c[~np.isnan(var_c)]
    bg = var_bg[~np.isnan(var_bg)]

    #create bins
    if k in _binning:
        bin_min, bin_max = _binning.get(k)
    else:
        max_sig, min_sig = sig.max(), sig.min()
        max_c, min_c = c.max(), c.min()
        max_bg, min_bg = bg.max(), bg.min()
        bin_max, bin_min = max(max_sig, max_c, max_bg), min(min_sig, min_c, min_bg)
    bins = np.linspace(bin_min, bin_max, 101)

    #histogram
    hist_sig, bins = np.histogram(sig, normed=True, bins=bins)
    hist_c, bins = np.histogram(c, normed=True, bins=bins)
    hist_bg, bins = np.histogram(bg, normed=True, bins=bins)

    histo_sig_collector.append(hist_sig)
    histo_c_collector.append(hist_c)
    histo_bg_collector.append(hist_bg)
    bin_collector.append(bins)

histo_sig_collector = np.asarray(histo_sig_collector)
histo_c_collector = np.asarray(histo_c_collector)
histo_bg_collector = np.asarray(histo_bg_collector)
bin_collector = np.asarray(bin_collector)

np.savetxt("histo_sig_collector.csv", histo_sig_collector, delimiter=',')
np.savetxt("histo_c_collector.csv", histo_c_collector, delimiter=',')
np.savetxt("histo_bg_collector.csv", histo_bg_collector, delimiter=',')
np.savetxt("bin_collector.csv", bin_collector, delimiter=',')
