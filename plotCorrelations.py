from matplotlib import pyplot as plt
from matplotlib import cm as cm
from pandas import DataFrame
import h5py


f = h5py.File('../' + 'gjj_Variables.hdf5', 'r')

high = f['high_input'][0:10000000]
highdf = DataFrame(high[:, 0, :])

fig = plt.figure()
ax1 = fig.add_subplot(111)
cmap = cm.get_cmap('jet', 30)
cax = ax1.imshow(highdf.corr(), cmap=cmap)
plt.title('Expert Level Feature Correlation')
high_var = ['jet_pt', 'jet_eta',
'track_2_d0_significance', 'track_3_d0_significance',
'track_2_z0_significance', 'track_3_z0_significance',
'n_tracks_over_d0_threshold', 'jet_prob', 'jet_width_eta', 'jet_width_phi',
'vertex_significance', 'n_secondary_vertices', 'n_secondary_vertex_tracks',
'delta_r_vertex', 'vertex_mass', 'vertex_energy_fraction']
plt.xticks(range(len(high_var)), high_var, fontsize=10, rotation='vertical')
plt.yticks(range(len(high_var)), high_var, fontsize=10)
fig.colorbar(cax)
plt.tight_layout()
plt.savefig('Expert_Correlations.png')
