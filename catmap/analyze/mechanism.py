from analysis_base import *
import math as np
from catmap.functions import convert_formation_energies

class MechanismAnalysis(MechanismPlot,ReactionModelWrapper,MapPlot):
    def __init__(self,reaction_model=None):
        self._rxm = reaction_model
        defaults = {'pressure_correction':True,
                    'min_pressure':1e-12,
                    'energy_type':'free_energy',
                    'include_labels':False,
                    'subplots_adjust_kwargs':{},
                    'kwarg_dict':{}}
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
        if not surfaces:
            surfaces = self.surface_names
        if not self.surface_colors:
            self.surface_colors = get_colors(max(len(surfaces),len(mechanisms)))
       
        self.kwarg_list = []
        for key in self.rxn_mechanisms.keys():
            self.kwarg_list.append(self.kwarg_dict.get(key,{}))

        for n,mech in enumerate(mechanisms):
            for i, surf in enumerate(surfaces):
                xy = self.descriptor_dict[surf]
                if '-' not in xy:
                    self.thermodynamics.current_state = None #force recalculation

                    if self.energy_type == 'free_energy':
                        energy_dict = self.scaler.get_free_energies(xy)

                    elif self.energy_type == 'potential_energy':
                        energy_dict = self.scaler.get_free_energies(xy)
                        energy_dict.update(self.scaler.get_electronic_energies(xy))

                    elif self.energy_type == 'interacting_energy':
                        if not self.interacting_energy_map:
                            raise ValueError('No interacting energy map found.')
                        G_dict = {}
                        G_labels = self.output_labels['interacting_energy']
                        valid = False
                        for pt, energies in self.interacting_energy_map:
                            if pt == xy:
                                valid = True
                                for ads,E in zip(G_labels, energies):
                                    G_dict[ads] = E
                        if valid == False:
                            raise UserWarning('No coverages found for '+xy+' in map')

                        if not G_dict:
                            raise ValueError('No energies found for point: ', xy)

                        energy_dict = self.scaler.get_free_energies(xy) #get gas G's
                        energy_dict.update(G_dict)

                    if self.pressure_correction == True:
                        for key in energy_dict:
                            if key.endswith('_g'):
                                P = self.gas_pressures[self.gas_names.index(key)]
                                energy_dict[key] += self._kB*self.temperature*np.log(P)
                   
                    if self.coverage_correction == True:
                        if not self.coverage_map:
                            raise ValueError('No coverage map found.')
                        cvg_labels = self.output_labels['interacting_energy']
                        valid = False
                        for pt, cvgs in self.coverage_map:
                            if pt == xy:
                                valid = True
                                for ads,cvg in zip(cvg_labels, cvgs):
                                    print ads, cvg
                                    energy_dict[ads] += self._kB*self.temperature*np.log(
                                                                                     cvg)
                        if valid == False:
                            raise UserWarning('No coverages found for '+xy+' in map')
                    
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

                            species = self.elementary_rxns[step-1][0]
                            L = self.label_maker(species)
                            self.labels.append(L)

                        else:
                            nrg = p[0]
                            bar = p[1]

                            species = self.elementary_rxns[step-1][-1]
                            L = self.label_maker(species)
                            self.labels.append(L)

                        self.energies.append(nrg)
                        self.barriers.append(bar)
                    
                    if labels and self.include_labels:
                        self.labels = labels
                    elif self.labels and self.include_labels:
                        pass
                    else:
                        self.labels = []

                    for e, e_a,rxn in zip(self.energies[1:],self.barriers,mech):
                        if rxn < 0:
                            reverse = True
                            rxn = rxn*-1
                        else:
                            reverse = False
                    self.data_dict[self.rxn_mechanisms.keys()[n]] = [self.energies,
                            self.barriers]
                    
                    kwargs = self.kwarg_list[n]
                    for key in kwargs:
                        setattr(self,key,kwargs[key])

                    self.draw(ax)
        
        if self.energy_type == 'free_energy':
            ax.set_ylabel('$\Delta G$ [eV]')
        elif self.energy_type == 'potential_energy':
            ax.set_ylabel('$\Delta E$ [eV]')
        if self.energy_type == 'interacting_energy':
            ax.set_ylabel('$\Delta G_{interacting}$ [eV]')
        fig.subplots_adjust(**self.subplots_adjust_kwargs)
        MapPlot.save(self,fig,
                save=save,default_name=self.model_name+'_pathway.pdf')
        self._fig = fig
        self._ax = ax
        return fig

    def label_maker(self,species):
        species = [s for s in species if self.species_definitions[s].get('type',None) not in 'site']
        new_species = []
        for sp in species:
            name,site = sp.split('_')
            for kj in range(0,9):
                name = name.replace(str(kj), '$_{'+str(kj)+'}$')
            if site == 'g':
                new_species.append(name+'(g)')
            else:
                new_species.append(name+'*')
        species_set = list(set(new_species))
        if species_set != new_species:
            species = [str(new_species.count(sp))+sp for sp in species_set]
        else:
            species = new_species

        L = '+'.join(species)
        return ' '+L
