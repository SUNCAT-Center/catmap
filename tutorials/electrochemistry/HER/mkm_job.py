from catmap import ReactionModel
import sys
from string import Template

voltage = sys.argv[1]

mkm_template = Template(open('HER_template.mkm').read())
mkm_text = mkm_template.substitute(voltage=voltage)
with open('HER.mkm','w') as f:
	f.write(mkm_text)
model = ReactionModel(setup_file = 'HER.mkm')
model.output_variables+=['production_rate', 'free_energy', 'selectivity']
model.run()

from catmap import analyze
vm = analyze.VectorMap(model)

vm.plot_variable = 'rate'
vm.log_scale = True
vm.min = 1e-10
vm.max = 1e6
fig = vm.plot(save=False)
fig.suptitle(str(voltage) + "V vs. RHE")
fig.savefig('rate' + voltage + '.png')

vm.plot_variable = 'production_rate'
vm.log_scale = True
vm.min = 1e-10
vm.max = 1e6
fig = vm.plot(save=False)
fig.suptitle(str(voltage) + "V vs. RHE")
fig.savefig('production_rate' + voltage + '.png')

vm.plot_variable = 'coverage'
vm.min = 0
vm.max = 1
vm.log_scale = False
fig = vm.plot(save=False)
fig.suptitle(str(voltage) + "V vs. RHE")
fig.savefig('coverage' + voltage + '.png')
