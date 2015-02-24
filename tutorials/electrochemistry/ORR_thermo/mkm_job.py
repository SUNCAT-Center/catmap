from catmap import ReactionModel
from string import Template

mkm_file = 'ORR.mkm'
model = ReactionModel(setup_file = mkm_file)
model.output_variables+=['production_rate', 'free_energy', 'selectivity']
model.run()

from catmap import analyze

vm = analyze.VectorMap(model)
vm.plot_variable = 'rate'
vm.log_scale = True
vm.min = 1e-6
vm.max = 1e4
fig = vm.plot(save='rate.png')

vm = analyze.VectorMap(model)
vm.plot_variable = 'production_rate'
vm.log_scale = True
vm.min = 1e-3
vm.max = 1e4
fig = vm.plot(save='prodrate.png')

vm = analyze.VectorMap(model)
vm.plot_variable = 'coverage'
vm.min = 0
vm.max = 1
vm.log_scale = False
fig = vm.plot(save='coverage.png')

# vm.model_summary()