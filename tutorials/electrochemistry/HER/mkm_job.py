from catmap import ReactionModel
import sys
from string import Template

model = ReactionModel(setup_file = 'HER.mkm')
model.output_variables+=['production_rate', 'free_energy', 'selectivity', 'directional_rates']
model.run()

from catmap import analyze
vm = analyze.VectorMap(model)

vm.plot_variable = 'rate'
vm.log_scale = True
vm.min = 1e-10
vm.max = 1e6
fig = vm.plot(save='rate.png')

vm.plot_variable = 'production_rate'
vm.log_scale = True
vm.min = 1e-10
vm.max = 1e6
fig = vm.plot(save='production_rate.png')

vm.plot_variable = 'coverage'
vm.min = 0
vm.max = 1
vm.log_scale = False
fig = vm.plot(save='coverage.png')

vm.plot_variable = 'directional_rates'
vm.log_scale = True
vm.min = 1e-20
vm.max = 1e-5
vm.unique_only = False  # Save all the figures (not just the unique one)
fig = vm.plot(save='directional_rates.png')