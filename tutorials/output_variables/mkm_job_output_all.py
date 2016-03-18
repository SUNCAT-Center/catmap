"""This mkm_job script generates the same results as the tutorialsoutput_variables/mkm_job.py script
however it demonstrates that values of output_variables do not need to be all hard-coded

FIXME: At this point we still need to hard-code two-dimensional variables, since this is difficult to glean from code inspection only.
FIXME: Output variables 'turnover_frequency' and 'interaction_matrix' currently crash, we should fix them or remove them from the parser.

"""
from catmap import ReactionModel
from catmap import analyze

mkm_file = 'CO_oxidation.mkm'
model = ReactionModel(setup_file=mkm_file)
model.output_variables = ['all']  # <==== Here is where the 'all' magic happens
two_dim_variables = ['rate_control', 'selectivity_control',
                     'rxn_order', 'frequency', 'interaction_matrix']
# !!! However we still need to hard-code those variables that are supposed to go through MatrixMap instead of VectorMap
model.run()

vm = analyze.VectorMap(model)
mm = analyze.MatrixMap(model)

def set_defaults(vm, mm):
    vm.log_scale = True
    vm.min = 1e-25
    vm.max = 1e3
    vm.threshold = 1e-25
    vm.unique_only = False
    mm.log_scale = False
    mm.min = -2
    mm.max = 2
    mm.unique_only = False

for output_variable in model.output_variables:
    m = mm if output_variable in two_dim_variables else vm
    plot_type = 'MatrixMap' if output_variable in two_dim_variables else 'VectorMap'

    if output_variable is None:
        pass
    # put other adjustments here
    # elif output_variable == '...':
    #     m.<attr> = '...'
    elif output_variable == 'interaction_matrix':
        continue  # only works for model with adsorbate-adsorbate interaction
    else:  # set default values
        # not necessarily the right choice of parameters for all
        # output_variables
        set_defaults(vm, mm)
    print(
        "Plotting {output_variable}, plot type: {plot_type}".format(**locals()))
    m.plot_variable = output_variable
    fig = m.plot(save='ALL_{output_variable}.pdf'.format(**locals()))
    fig.clf()
