from analysis_base import *

class ScalingAnalysis(ScalingPlot,ReactionModelWrapper):
    """
    .. todo:: __doc__
    """
    def __init__(self,reaction_model):
        self._rxm = reaction_model
        self.scaling_mode = 'linear'
        ScalingPlot.__init__(self,self.descriptor_names,
                self.descriptor_dict,self.surface_names,
                self.parameter_dict,
                self.plotter_scaling_function,
                self.plotter_x_axis_function)

    def plotter_scaling_function(self,descriptors,**kwargs):
        """
        .. todo:: __doc__
        """
        descriptors = list(descriptors)
        return self.scaler.get_electronic_energies(descriptors,**kwargs)

    def plotter_x_axis_function(self,descriptors,**kwargs):
        """
        .. todo:: __doc__
        """
        if self.scaling_mode == 'linear':
            x_dict = {}
            labels = {}
            if self.coefficient_matrix is None:
                C = self.scaler.get_coefficient_matrix()
                self.coefficient_matrix = C
            else:
                C = self.coefficient_matrix
            for ads,c in zip(
                    self.adsorbate_names+self.transition_state_names,C):
                x_dict[ads] = sum(
                        [c[i]*a for i,a in enumerate(descriptors)]) + c[-1]
                lab = ''
                for i,a in enumerate(self.descriptor_names):
                    c_print =  round(c[i],self.descriptor_decimal_precision)
                    if c_print > 0:
                        lab += str(c_print)+'*$E_{'+str(a)+'}$+'
                const = round(c[-1],self.descriptor_decimal_precision)
                if const > 0:
                    lab += '+' + str(const)
                elif const < 0:
                    lab += '-' + str(-const)
                lab = lab.replace('++','+')
                lab = lab.replace('+-','-')
                labels[ads] = lab
            return x_dict,labels

    def get_error(self):
        """
        .. todo:: __doc__
        """
        if not self.scaling_error:
            self.plot(save=False)
        return self.scaling_error

