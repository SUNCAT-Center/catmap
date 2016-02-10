from catmap import ReactionModel

mkm_file = 'CO_oxidation.mkm'
model = ReactionModel(setup_file=mkm_file)
model.output_variables += ['production_rate']#,'rate_control']
model.run()

#model = ReactionModel(setup_file=mkm_file.replace('mkm','log'))

from catmap import analyze
vm = analyze.VectorMap(model)
vm.plot_variable = 'production_rate' #tell the model which output to plot
vm.log_scale = True #rates should be plotted on a log-scale
vm.min = 1e-25 #minimum rate to plot
vm.max = 1e3 #maximum rate to plot

vm.descriptor_labels = ['O reactivity [eV]', 'CO reactivity [eV]', ]
vm.threshold = 1e-25 #anything below this is considered to be 0
vm.subplots_adjust_kwargs = {'left':0.2,'right':0.8,'bottom':0.15}
vm.plot(save='pretty_production_rate.pdf')

#mm = analyze.MatrixMap(model)
#mm.plot_variable = 'rate_control'
#mm.log_scale = False
#mm.min = -2
#mm.max = 2
#mm.plot(save='rate_control.pdf')
