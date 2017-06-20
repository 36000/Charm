'''
high level variables (16):

[jet_pt, jet_eta,
track_2_d0_significance, track_3_d0_significance,
track_2_z0_significance, track_3_z0_significance,
n_tracks_over_d0_threshold, jet_prob, jet_width_eta, jet_width_phi,
vertex_significance, n_secondary_vertices, n_secondary_vertex_tracks,
delta_r_vertex, vertex_mass, vertex_energy_fraction]

flavor (y):
signal --> y == 5
'''

_binning = {0:[0,300], 1:[-3,3], 2:[0,2.5], 3:[0,5],
4:[0,5], 5:[0,5], 6:[0,10], 7:[0,0.04], 
8:[0,0.4], 9:[0,0.4], 10:[0,5], 11:[0,10], 
12:[0,10], 13:[0,7], 14:[0,25], 15:[0,5]}

high_var = ['jet_pt', 'jet_eta',
'track_2_d0_significance', 'track_3_d0_significance',
'track_2_z0_significance', 'track_3_z0_significance',
'n_tracks_over_d0_threshold', 'jet_prob', 'jet_width_eta', 'jet_width_phi',
'vertex_significance', 'n_secondary_vertices', 'n_secondary_vertex_tracks',
'delta_r_vertex', 'vertex_mass', 'vertex_energy_fraction']

import numpy as np
import h5py
import matplotlib.pyplot as plt

#filepath = '/phys/groups/tev/scratch1/users/jk232/'
#f = h5py.File(filepath+'gjj_Variables.hdf5', 'r')

f = h5py.File('../'+'gjj_Variables.hdf5', 'r')

high = f['high_input'][0:10000000]
y = f['y_input'][0:10000000]

high_sig_collect = high[y[:,2].astype(bool), :,:]
high_c_collect = high[y[:,1].astype(bool), :,:]
high_bg_collect = high[y[:,0].astype(bool),:,:]

# same as Rev1 below
histo_sig_collector = []
histo_bg_collector = []
histo_c_collector = []
bin_collector = []

#histogram for each variable k
for k in range(high_sig_collect.shape[2]):

    var_sig = high_sig_collect[:,:,k]
    var_c = high_c_collect[:,:,k]
    var_bg = high_bg_collect[:,:,k]

    #remove nan for plotting
    sig = var_sig[~np.isnan(var_sig)]
    bg = var_bg[~np.isnan(var_bg)]
    c = var_c[~np.isnan(var_c)]

    #create bins
    if k in _binning:
        bin_min, bin_max = _binning.get(k)
    else:
        max_sig, min_sig = sig.max(), sig.min()
        max_bg, min_bg = bg.max(), bg.min()
        max_c, min_c = c.max(), c.min()
        bin_max, bin_min = max(max_sig, max_bg, max_c), min(min_sig, min_bg, min_c)
    if k == 6 or  k == 11 or  k == 12:
        bins = np.linspace(bin_min, bin_max, 11)
    else:
        bins = np.linspace(bin_min, bin_max, 101)

    #histogram
    hist_sig, bins = np.histogram(sig, normed=True, bins=bins)
    hist_c, bins = np.histogram(c, normed=True, bins=bins)
    hist_bg, bins = np.histogram(bg, normed=True, bins=bins)
    
    pltBg = plt.plot(bins[:-1], hist_bg, drawstyle='steps-post', color='0.2', label='bg')
    pltC = plt.plot(bins[:-1], hist_c, drawstyle='steps-post', color='blue', label='charm')
    pltSig = plt.plot(bins[:-1], hist_sig, drawstyle='steps-post', color='red', label='bottom')
    if k != 1:
        plt.yscale('log')
    plt.title(high_var[k])
    plt.legend(loc='upper right')
    plt.savefig('h'+str(k)+'.png')
    plt.clf()

histo_sig_collector = np.asarray(histo_sig_collector)
histo_c_collector = np.asarray(histo_c_collector)
histo_bg_collector = np.asarray(histo_bg_collector)
bin_collector = np.asarray(bin_collector)

np.savetxt("histo_sig_collector.csv", histo_sig_collector, delimiter=',')
np.savetxt("histo_c_collector.csv", histo_c_collector, delimiter=',')
np.savetxt("histo_bg_collector.csv", histo_bg_collector, delimiter=',')
np.savetxt("bin_collector.csv", bin_collector, delimiter=',')
