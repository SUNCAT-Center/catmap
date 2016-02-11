from catmap import ReactionModel

mkm_file = 'CO_oxidation.mkm'
model = ReactionModel(setup_file=mkm_file)
scaler_variables = ['rxn_parameter','frequency','electronic_energy',
    'free_energy','zero_point_energy','enthalpy','entropy']
# Warning: adding 'rate_control' or 'selectivity_control' to output_variables increases calculation time significantly
solver_variables = ['production_rate','rate_control','coverage','rate',
    'consumption_rate','selectivity','rxn_direction',
    'rate_constant','equilibrium_constant','carbon_selectivity',
    'selectivity_control','rxn_order','apparent_activation_energy',
    'interacting_energy','directional_rates','forward_rate','reverse_rate',
    'forward_rate_constant','reverse_rate_constant']
two_dim_variables = ['rate_control', 'selectivity_control', 'rxn_order', 'frequency']
one_dim_variables = list(set(scaler_variables + solver_variables) - set(two_dim_variables))
model.output_variables += solver_variables + scaler_variables

model.run()

from catmap import analyze
vm = analyze.VectorMap(model)
vm.log_scale = True  # not necessarily the right choice of parameters for all output_variables
vm.min = 1e-25
vm.max = 1e3
vm.threshold = 1e-25
vm.unique_only = False
for out_var in one_dim_variables:
    vm.plot_variable = out_var
    fig = vm.plot(save=out_var + '.pdf')
    fig.clf()

mm = analyze.MatrixMap(model)
mm.log_scale = False
mm.min = -2
mm.max = 2
mm.unique_only = False
for out_var in two_dim_variables:
    mm.plot_variable = out_var
    fig = mm.plot(save=out_var + '.pdf')
    fig.clf()
