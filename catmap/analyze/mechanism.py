from analysis_base import *
import numpy as np
from math import log
from catmap.functions import convert_formation_energies

class MechanismAnalysis(MechanismPlot,ReactionModelWrapper,MapPlot):
    """
      A simple tool for the generation of potential energy diagrams
      from a reaction network.
    """
    def __init__(self,reaction_model=None):
        """
        Class for generating potential energy diagrams.

        :param reaction_model: The ReactionModel object to load.
        """
        self._rxm = reaction_model
        defaults = {'pressure_correction':True,
                    'min_pressure':1e-12,
                    'energy_type':'free_energy',
                    'include_labels':False,
                    'subplots_adjust_kwargs':{},
                    'line_offset':0,
                    'kwarg_dict':{}}
        self._rxm.update(defaults)
        self.data_dict = {}
        MechanismPlot.__init__(self,[0])

    def plot(self,ax=None,plot_variants=None,mechanisms=None,
            labels=None,save=True):
        """
        Generates the potential energy diagram plot

        :param ax: Matplotlib Axes object to optionally plot into

        :param plot_variants: Which PEDs to plot. Defaults to all surfaces
                              or all applied voltages
        :param plot_variants: list of voltages (if electrochemical) or
                              descriptor tuples to plot

        :param mechanisms: Which reaction pathways to plot.  Each integer
                           corresponds to an elementary step. Elementary
                           steps are indexed in the order that they are
                           input with 1 being the first index. Negative
                           integers are used to designate reverse reactions.
                           Read in from model.rxn_mechanisms by default
        :type mechanisms: {string:[int]}

        :param labels: Labels for each state.  Can be generated automatically
        :type labels: [string]

        :param save: Whether to write plot to file
        :type save: bool
        """
        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            fig = ax.get_figure()
        if not mechanisms:
            mechanisms = self.rxn_mechanisms.values()
        if not plot_variants and self.descriptor_dict:
            plot_variants = self.surface_names
        elif not plot_variants and 'voltage' in self.descriptor_names:
            voltage_idx = self.descriptor_names.index('voltage')
            v_min, v_max = self.descriptor_ranges[voltage_idx]
            plot_variants = np.linspace(v_min, v_max, 5)
        if not self.plot_variant_colors:
            self.plot_variant_colors = get_colors(max(len(plot_variants),len(mechanisms)))

        self.kwarg_list = []
        for key in self.rxn_mechanisms.keys():
            self.kwarg_list.append(self.kwarg_dict.get(key,{}))

        for n,mech in enumerate(mechanisms):
            for i, variant in enumerate(plot_variants):
                if self.descriptor_dict:
                    xy = self.descriptor_dict[variant]
                elif 'voltage' in self.descriptor_names:
                    voltage_idx = self.descriptor_names.index('voltage')
                    xy = [0, 0]
                    xy[voltage_idx] = variant
                    xy[1-voltage_idx] = self.descriptor_ranges[1-voltage_idx][0]
                else:
                    xy = variant
                if '-' not in xy:
                    self.thermodynamics.current_state = None #force recalculation
                    self._rxm._descriptors = xy

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
                        xyo = self.nearest_mapped_point(self.interacting_energy_map,xy)
                        if xyo != xy:
                            print('Warning: No output at: '+str(xy)+'. Using output from: '+str(xyo))
                        xy = xyo
                        self._rxm._descriptors = xyo

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
                                energy_dict[key] += self._kB*self.temperature*log(P)
                   
                    if self.coverage_correction == True:
                        if not self.coverage_map:
                            raise UserWarning('No coverage map found.')
                        cvg_labels = self.output_labels['interacting_energy']
                        valid = False
                        for pt, cvgs in self.coverage_map:
                            if pt == xy:
                                valid = True
                                for ads,cvg in zip(cvg_labels, cvgs):
                                    energy_dict[ads] += self._kB*self.temperature*log(
                                                                            float(cvg))
                        if valid == False:
                            raise UserWarning('No coverages found for '+str(xy)+' in map')
                    
                    params = self.adsorption_to_reaction_energies(energy_dict)
                    self.energies = [0]
                    self.barriers = []
                    self.labels = ['']
                    if len(plot_variants) > 1:
                        self.energy_line_args['color'] = \
                                self.barrier_line_args['color'] = \
                                self.plot_variant_colors[i]
                    else:
                        self.energy_line_args['color'] = \
                                self.barrier_line_args['color'] = \
                                self.plot_variant_colors[n]
                    for step in mech:

                        if str(step).startswith('half'):
                            step = int(step.replace('half',''))
                            split = True
                        else:
                            split = False

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
                            if split:
                                L = L.strip()
                                if L.startswith('2'):
                                    L = L[1:]
                                else:
                                    L = r'$\frac{1}{2}$'+L
                                L = ' '+L #add padding back.
                            self.labels.append(L)

                        if split == False:
                            self.energies.append(nrg)
                            self.barriers.append(bar)
                        elif split == True:
                            self.energies.append(0.5*nrg)
                            self.barriers.append(0) #split steps cannot have barriers.

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
                    
                    self.initial_energy += self.line_offset
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
        """
        .. todo:: __doc__
        """
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
