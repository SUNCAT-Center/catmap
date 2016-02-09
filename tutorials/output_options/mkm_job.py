from catmap import ReactionModel

mkm_file = 'CO_oxidation.mkm'
model = ReactionModel(setup_file=mkm_file)
model.output_variables += ['production_rate','rate_control','coverage','rate',
    'consumption_rate','selectivity','rxn_direction','selectivity',
    'rate_constant','equilibrium_constant','carbon_selectivity',
    'selectivity_control','rxn_order','apparent_activation_energy',
    'interacting_energy','directional_rates','forward_rate','reverse_rate',
    'forward_rate_constant','reverse_rate_constant',
    'rxn_parameter','frequency','electronic_energy',
    'free_energy','zero_point_energy','enthalpy','entropy']
model.run()

#model = ReactionModel(setup_file=mkm_file.replace('mkm','log'))

from catmap import analyze
vm = analyze.VectorMap(model)
vm.plot_variable = 'production_rate' #tell the model which output to plot
vm.log_scale = True #rates should be plotted on a log-scale
vm.min = 1e-25 #minimum rate to plot
vm.max = 1e3 #maximum rate to plot
vm.unique_only = False

vm.descriptor_labels = ['CO reactivity [eV]', 'O reactivity [eV]']
vm.threshold = 1e-25 #anything below this is considered to be 0
vm.subplots_adjust_kwargs = {'left':0.2,'right':0.8,'bottom':0.15}
vm.plot(save='pretty_production_rate.pdf')

for out_var in model.output_variables:
    if out_var in ['rate_control', 'production_rate', 'selectivity_control', 'rxn_order', 'frequency']:
        continue
    print(out_var)
    vm.plot_variable = out_var
    vm.plot(save=out_var + '.pdf')

# Warning: adding 'rate_control' to output_variables increases calculation time significantly
mm = analyze.MatrixMap(model)
for out_var in ['rate_control', 'selectivity_control', 'rxn_order', 'frequency']:
    print(out_var)
    mm.plot_variable = out_var
    mm.log_scale = False
    mm.min = -2
    mm.max = 2
    mm.plot(save=out_var + '.pdf')

