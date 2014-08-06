import pickle
from catmap.model import ReactionModel as RM
import numpy as np
import catmap
from glob import glob
import os
from matplotlib import pyplot as plt
import sys

"""
The first argument should be the name of the directory that contains all the various
subdirectories of catmap calculations at various voltages.  It should be the same
argument as was passed to make_plots.sh
"""

def get_data(which_data, minval=1e-5, base_glob='test11/*', rxn_index=1, log_scale=True):
	# Pt(111) OH binding energy
	Pt_x = 0.75

	dirs = glob(base_glob)
	voltages = []
	out_datas = []

	for base_dir in dirs:
		path = os.path.join(base_dir, 'ORR.pkl')
		voltage = base_dir.split('voltage')[-1]
		with open(path) as f:
			data = pickle.load(f)
			prod_rate_map = data[which_data]

			# from vectormap to make 1d data from 2d data
			pts,datas = zip(*prod_rate_map)
			pts = [[first] for first, second in pts]
			cols = zip(*datas)
			xy, rates = pts,cols
			mapp = zip(pts,datas)

			x = zip(*xy)[0]
			descriptor_ranges = [[min(x),max(x)]]
			eff_res = [len(xy), 1]
			maparray = RM.map_to_array(mapp,descriptor_ranges,eff_res,
			        log_interpolate=log_scale,minval=minval)

			x_range = descriptor_ranges[0]
			plot_in = [np.linspace(*x_range+eff_res),maparray[:,rxn_index]]
			y_sp = catmap.spline(plot_in[0], plot_in[1], k=1)
			x, y = Pt_x, y_sp(Pt_x)
			print base_dir, voltage, x, y
			voltages.append(voltage)
			out_datas.append(y)
	return voltages, out_datas

plots_directory = sys.argv[1]

# H2O production rate plots
if 1:
	voltages, H2O_rates = get_data('production_rate_map', minval=1e-5, base_glob=plots_directory + '/*', rxn_index=1, log_scale=True)

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(voltages, H2O_rates)
	ax.invert_yaxis()
	ax.set_xlabel('Potential (V vs RHE)')
	ax.set_ylabel(r'Production rate of H2O_g (s$^{-1}$)')
	fig.savefig('H2O_production.png')

# Coverage plots
if 1:
	voltages, OHA_coverage = get_data('coverage_map', minval=1e-5, base_glob=plots_directory + '/*', rxn_index=3, log_scale=True)
	voltages, OA_coverage = get_data('coverage_map', minval=1e-5, base_glob=plots_directory + '/*', rxn_index=5, log_scale=True)
	empty_coverage = [1 - OHA_coverage[i] - OA_coverage[i] for i in range(len(voltages))]

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax1, = ax.plot(voltages, OHA_coverage, 'y-')
	ax2, = ax.plot(voltages, OA_coverage, 'b-')
	ax3, = ax.plot(voltages, empty_coverage, 'r-')
	labels = ['OH*', 'O*', '1-OH*-O*']
	ax.set_ylim([-0.02, 1.02])
	ax.legend([ax1, ax2, ax3], labels, loc='best')
	ax.set_xlabel('Potential (V vs RHE)')
	ax.set_ylabel(r'Coverage')
	fig.savefig('coverage.png')