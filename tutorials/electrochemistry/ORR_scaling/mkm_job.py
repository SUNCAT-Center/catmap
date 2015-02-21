from catmap import ReactionModel
from string import Template


# do calculation
mkm_file = 'ORR.mkm'
model = ReactionModel(setup_file = mkm_file)
model.output_variables+=['production_rate', 'free_energy', 'selectivity']
model.run()

from catmap import analyze

vm = analyze.VectorMap(model)
vm.plot_variable = 'rate'
vm.log_scale = True
vm.min = 1e-10
vm.max = 1e5
vm.unique_only = False
fig = vm.plot(save="rate.png")

vm = analyze.VectorMap(model)
vm.plot_variable = 'production_rate'
vm.log_scale = True
vm.min = 1e-5
vm.max = 1e5
fig = vm.plot(save="prod_rate.png")

vm = analyze.VectorMap(model)
vm.plot_variable = 'coverage'
vm.min = 0
vm.max = 1
vm.log_scale = False
vm.unique_only = False
fig = vm.plot(save='coverage.png')

sa = analyze.ScalingAnalysis(model)
sa.plot(save='scaling.pdf')

vm.model_summary()