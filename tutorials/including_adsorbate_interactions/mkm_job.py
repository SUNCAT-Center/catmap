import pylab as plt
from catmap import ReactionModel

mkm_file = 'CO_oxidation.mkm'
model = ReactionModel(setup_file=mkm_file)
model.output_variables += ['production_rate','interaction_matrix','free_energy']
model.create_standalone = True
model.run()

from catmap import analyze
vm = analyze.VectorMap(model)
vm.plot_variable = 'production_rate' #tell the model which output to plot
vm.log_scale = True #rates should be plotted on a log-scale
vm.min = 1e-25 #minimum rate to plot
vm.max = 1e8 #maximum rate to plot
vm.threshold = 1e-25 #anything below this is considered to be 0
vm.subplots_adjust_kwargs = {'left':0.2,'right':0.8,'bottom':0.15}
vm.plot(save='production_rate.pdf')

vm.plot_variable = 'coverage' #tell the model which output to plot
vm.log_scale = False #coverage should not be plotted on a log-scale
vm.min = 0 #minimum coverage
vm.max = 1 #maximum coverage
vm.subplots_adjust_kwargs = {'left':0.1,'right':0.85,'wspace':0.6,'bottom':0.15}
vm.plot(save='coverage.pdf')
