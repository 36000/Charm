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


fig = plt.figure(0, figsize=(35,35))
fig.suptitle('ROC curves for jet level variables', fontsize=80, weight='roman')
fig.subplots_adjust(hspace=0.4, wspace=0.4)

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
    
    #normal plots
    plt.figure(1)
    pltBg = plt.plot(bins[:-1], hist_bg, drawstyle='steps-post', color='0.2', label='bg')
    pltC = plt.plot(bins[:-1], hist_c, drawstyle='steps-post', color='blue', label='charm')
    pltSig = plt.plot(bins[:-1], hist_sig, drawstyle='steps-post', color='red', label='bottom')
    if k != 1:
        plt.yscale('log')
    plt.title(high_var[k])
    plt.legend(loc='upper right')
    plt.savefig('h'+str(k)+'.png')
    plt.clf()
    
    #ROC
    tpr = []
    fpr = []
    tprc = []
    fprc = []

    sigFreqSum = np.sum(hist_sig)
    cFreqSum = np.sum(hist_c)
    bgFreqSum = np.sum(hist_bg)
    
    for j in range(hist_sig.shape[0]): # looping each histo_sig
        if high_var[k] in ['jet_prob', 'delta_r_vertex']:
            TP = np.sum(hist_sig[0:j+1])
            FN = np.sum(hist_bg[0:j+1])
            tpr.append(TP / float(sigFreqSum))
            fpr.append(FN / float(bgFreqSum))
        else:    # flip symmetry along dash line in ROC
            TP = np.sum(hist_sig[j:])
            FN = np.sum(hist_bg[j:])
            tpr.append(TP / float(sigFreqSum))
            fpr.append(FN / float(bgFreqSum))
            
    for j in range(hist_c.shape[0]): # looping each histo_sig
        if high_var[k] in ['jet_prob', 'delta_r_vertex']:
            TP = np.sum(hist_c[0:j+1])
            FN = np.sum(hist_bg[0:j+1])
            tprc.append(TP / float(cFreqSum))
            fprc.append(FN / float(bgFreqSum))
        else:    # flip symmetry along dash line in ROC
            TP = np.sum(hist_c[j:])
            FN = np.sum(hist_bg[j:])
            tprc.append(TP / float(cFreqSum))
            fprc.append(FN / float(bgFreqSum))

    if high_var[k] in ['jet_prob', 'delta_r_vertex']:
        tpr = [0.0] + tpr
        fpr = [0.0] + fpr
        tprc = [0.0] + tprc
        fprc = [0.0] + fprc

    tpr, fpr, tprc, fprc = np.asarray(tpr), np.asarray(fpr), np.asarray(tprc), np.asarray(fprc)
    
    #AUC
    ufpr = np.unique(fpr)
    utpr = []
    for i in range(len(ufpr)):
        ind = np.where(fpr == ufpr[i])
        utpr.append(np.max(tpr[ind]))

    utpr = np.asarray(utpr)
    
    AUCb = np.trapz(utpr, x=ufpr)
    st_AUCb = "AUC = " + str('%.2f' % round(AUCb, 2))
    
    ufprc = np.unique(fprc)
    utprc = []
    for i in range(len(ufprc)):
        ind = np.where(fprc == ufprc[i])
        utprc.append(np.max(tprc[ind]))

    utprc = np.asarray(utprc)

    AUCc = np.trapz(utprc, x=ufprc)
    st_AUCc = "AUC = " + str('%.2f' % round(AUCc, 2))
    
    plt.figure(0)
    ax = fig.add_subplot(4,4,k+1)
    ax.plot(fprc, tprc, linewidth=2, color='blue', lw=2)   #drawstyle='steps-post'
    ax.plot(fpr, tpr, linewidth=2, color='red', lw=2)   #drawstyle='steps-post'
    ax.plot([0,1],[0,1], color='0.2', lw=1.5)
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.text(0.55,0.15, st_AUCb, fontsize=23, weight=550)
    ax.text(0.55,0.05, st_AUCc, fontsize=23, weight=550)
    ax.set_xlabel('false positive rate (fpr)', fontsize=27)
    ax.set_ylabel('true positive rate (tpr)', fontsize=27)
    ax.set_title(high_var[k], fontsize=30, weight=550)

fig.savefig('High_Level_PR_Rev3'+'.png')

