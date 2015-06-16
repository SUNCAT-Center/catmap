from catmap import ReactionModel
from catmap import analyze
import numpy as np
import cPickle as pickle

include_overbinding = False
include_rate_control = False
mkm_file = 'EtOH.mkm'
model = ReactionModel(setup_file=mkm_file)
model.create_standalone = True

if include_overbinding:
    overbinding = 0.25
    CO_energies = model.species_definitions['CO_s']['formation_energy']
    CO_energies = [E+overbinding for E in CO_energies]
    model.species_definitions['CO_s']['formation_energy'] = CO_energies

model.output_variables += ['production_rate','selectivity','zero_point_energy','enthalpy','entropy','free_energy']
model.output_variables += ['interaction_matrix','interacting_energy','equilibrium_constant']
if include_rate_control:
    model.output_variables += ['rate_control']

model.run()

vm = analyze.VectorMap(model)
vm.log_scale = True #rates should be plotted on a log-scale
vm.min = 1e-15 #minimum rate to plot
vm.max = 1 #maximum rate to plot
vm.threshold = 1e-25 #do not plot rates below this
vm.descriptor_labels = ['Temperature [K]','log$_{10}$(Pressure [bar])']
vm.plot_variable = 'production_rate'
vm.plot(save='production_rate.pdf')

vm.plot_variable = 'coverage'
vm.log_scale = False
vm.min = 0
vm.max = 1
vm.plot(save='coverage.pdf')

vm.plot_variable = 'selectivity'
vm.subplots_adjust_kwargs['wspace'] = 0.55
vm.plot(save='selectivity.pdf')

if include_rate_control:
    mm = analyze.MatrixMap(model)
    mm.plot_variable = 'rate_control'
    mm.log_scale = False
    mm.min = -2
    mm.max = 2
    mm.plot(save='rate_control.pdf')
