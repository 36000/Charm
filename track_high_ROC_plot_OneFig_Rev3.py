high_var = ['jet_pt', 'jet_eta',
'track_2_d0_significance', 'track_3_d0_significance',
'track_2_z0_significance', 'track_3_z0_significance',
'n_tracks_over_d0_threshold', 'jet_prob', 'jet_width_eta', 'jet_width_phi',
'vertex_significance', 'n_secondary_vertices', 'n_secondary_vertex_tracks',
'delta_r_vertex', 'vertex_mass', 'vertex_energy_fraction']


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps

# input files has dimensions [variables x frequencies]
histo_sig = np.genfromtxt("histo_sig_collector.csv", delimiter=',')
histo_bg = np.genfromtxt("histo_bg_collector.csv", delimiter=',')
bino = np.genfromtxt("bin_collector.csv", delimiter=',')

fig = plt.figure(figsize=(35,35))
fig.suptitle('ROC curves for jet level variables', fontsize=80, weight='roman')
fig.subplots_adjust(hspace=0.4, wspace=0.4)

for i in range(histo_sig.shape[0]): # looping variables

    tpr = []
    fpr = []

    sigFreqSum = np.sum(histo_sig[i,:])
    bgFreqSum = np.sum(histo_bg[i,:])

    for j in range(histo_sig.shape[1]): # looping each histo_sig
        if high_var[i] in ['jet_prob', 'delta_r_vertex']:
            TP = np.sum(histo_sig[i, 0:j+1])
            FN = np.sum(histo_bg[i, 0:j+1])
            tpr.append(TP / float(sigFreqSum))
            fpr.append(FN / float(bgFreqSum))
        else:    # flip symmetry along dash line in ROC
            TP = np.sum(histo_sig[i,j:])
            FN = np.sum(histo_bg[i,j:])
            tpr.append(TP / float(sigFreqSum))
            fpr.append(FN / float(bgFreqSum))

    if high_var[i] in ['jet_prob', 'delta_r_vertex']:
        tpr = [0.0] + tpr
        fpr = [0.0] + fpr

    tpr, fpr = np.asarray(tpr), np.asarray(fpr)

    #AUC
    ufpr = np.unique(fpr)
    utpr = []
    for k in range(len(ufpr)):
        ind = np.where(fpr == ufpr[k])
        utpr.append(np.max(tpr[ind]))

    utpr = np.asarray(utpr)

    AUC = np.trapz(utpr, x=ufpr)
    st_AUC = "AUC = " + str('%.2f' % round(AUC, 2))

    ax = fig.add_subplot(4,4,i+1)
    ax.plot(fpr, tpr, linewidth=2, color='blue', lw=2)   #drawstyle='steps-post'
    ax.plot([0,1],[0,1], 'r--', lw=1.5)
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.text(0.55,0.07, st_AUC, fontsize=23, weight=550)
    ax.set_xlabel('false positive rate (fpr)', fontsize=27)
    ax.set_ylabel('true positive rate (tpr)', fontsize=27)
    ax.set_title(high_var[i], fontsize=30, weight=550)

fig.savefig('High_Level_PR_Rev3'+'.png')
