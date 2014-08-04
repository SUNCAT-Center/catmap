from analysis_base import *
import math as np

class MechanismAnalysis(MechanismPlot,ReactionModelWrapper,MapPlot):
    def __init__(self,reaction_model=None):
        self._rxm = reaction_model
        defaults = {'pressure_correction':True,
                    'min_pressure':1e-12}
        self._rxm.update(defaults)
        self.data_dict = {}
        MechanismPlot.__init__(self,[0])

    def plot(self,ax=None,surfaces=None,mechanisms=None,
            labels=None,save=True):
        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            fig = None
        if not mechanisms:
            mechanisms = self.rxn_mechanisms.values()
            if not labels:
                labels = self.rxn_mechanisms.keys()
        if not surfaces:
            surfaces = self.surface_names
        if not self.surface_colors:
            self.surface_colors = get_colors(max(len(surfaces),len(mechanisms)))
        for n,mech in enumerate(mechanisms):
            for i, surf in enumerate(surfaces):
                xy = self.descriptor_dict[surf]
                if '-' not in xy:
                    self.thermodynamics.current_state = None #force recalculation
                    energy_dict = self.scaler.get_free_energies(xy)
                    params = self.adsorption_to_reaction_energies(energy_dict)

                    self.energies = [0]
                    self.barriers = []
                    self.labels = ['']
                    if len(surfaces) > 1:
                        self.energy_line_args['color'] = \
                                self.barrier_line_args['color'] = \
                                self.surface_colors[i]
                    else:
                        self.energy_line_args['color'] = \
                                self.barrier_line_args['color'] = \
                                self.surface_colors[n]
                    for step in mech:
                        if step < 0: 
                            reverse = True
                            step = step*-1
                        else: 
                            reverse = False
                        p = params[step-1]
                        if reverse == True:
                            nrg = -1*p[0]
                            bar = p[1] + nrg
                            self.labels.append('+'.join(self.elementary_rxns[step-1][0]))
                        else:
                            nrg = p[0]
                            bar = p[1]
                            self.labels.append('+'.join(self.elementary_rxns[step-1][-1]))
                        correction = 0
                        if self.pressure_correction == True:
                            IS_gasses = [self.gas_pressures[self.gas_names.index(s)]
                                    for s in self.elementary_rxns[step-1][0] 
                                    if s in self.gas_names]
                            FS_gasses = [self.gas_pressures[self.gas_names.index(s)]
                                    for s in self.elementary_rxns[step-1][-1] 
                                    if s in self.gas_names]
                            P_IS = max(sum(IS_gasses),self.min_pressure)
                            P_FS = max(sum(FS_gasses),self.min_pressure)
                            if IS_gasses:
                                correction -= float(
                                        self._kB*self.temperature*np.log(P_IS))
                                if bar != 0:
                                    bar = max(0,bar+correction)
                            if FS_gasses:
                                correction += float(
                                        self._kB*self.temperature*np.log(P_FS))
                            if reverse:
                                correction = -correction
                        self.energies.append(nrg+correction)
                        self.barriers.append(bar)

                    if labels:
                        self.labels = labels

                    for e, e_a,rxn in zip(self.energies[1:],self.barriers,mech):
                        if rxn < 0:
                            reverse = True
                            rxn = rxn*-1
                        else:
                            reverse = False
                    self.data_dict[self.rxn_mechanisms.keys()[n]] = [self.energies,
                            self.barriers]
                    self.draw(ax)
        MapPlot.save(self,fig,
                save=save,default_name=self.model_name+'_pathway.pdf')
        self._fig = fig
        self._ax = ax
        return fig
