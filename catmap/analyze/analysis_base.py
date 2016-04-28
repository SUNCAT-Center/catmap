import catmap
from catmap import ReactionModelWrapper
from catmap.model import ReactionModel as RM
from catmap import griddata
from copy import copy
try:
    from scipy.stats import norm
except:
    norm = None
from matplotlib.ticker import MaxNLocator
import os
import math
plt = catmap.plt
pickle= catmap.pickle
np = catmap.np
spline = catmap.spline
mtransforms = catmap.mtransforms

basic_colors = [[0,0,0],[0,0,1],[0.1,1,0.1],[1,0,0],[0,1,1],[1,0.5,0],[1,0.9,0],
                [1,0,1],[0,0.5,0.5],[0.5,0.25,0.15],[0.5,0.5,0.5]]
               #black,blue,green,red,cyan,orange,yellow,magenta,turquoise,brown,gray

def get_colors(n_colors):
    """
    Get n colors from basic_colors.

    :param n_colors: Number of colors
    :type n_colors: int
    """
    if n_colors <len(basic_colors):
        return basic_colors[0:n_colors]
    else:
        longlist= basic_colors*n_colors
        return longlist[0:n_colors]

def boltzmann_vector(energy_list,vector_list,temperature):
    """
    Create a vector which is a Boltzmann average of the vector_list weighted
    with energies in the energy_list.

    :param energy_list: List of energies
    :type energy_list: list

    :param vector_list: List of vectors
    :type energy_list: list

    :param temperature: Temperature
    :type energy_list: float
    """
    def boltzmann_avg(es,ns,T):
        """
        Calculate the Boltzmann average
        
        :param es: energies
        :type es: iterable

        :param ns:
        :type ns: iterable

        :param T: temperature
        :type T: float

        ..todo: description for ns
        """
        kB = 8.613e-5 #assuming energies are in eV and T is in K
        es = [e-min(es) for e in es] #normalize to minimum energy
        exp_sum = sum([np.exp(-e/(kB*T)) for e in es])
        exp_weighted = [n*np.exp(-e/(kB*T))/exp_sum for n,e in zip(ns,es)]
        Z = sum(exp_weighted)
        return Z
    vars = zip(*vector_list)
    boltz_vec = [boltzmann_avg(energy_list,v,temperature) for v in vars]
    return boltz_vec

class MapPlot:
    """
    Class for generating plots using a dictionary of default plotting attributes.

    The following attributes can be modified:

    :param resolution_enhancement: Resolution enhancement for interpolated maps
    :type resolution_enhancement: int

    :param min: Minimum
    :type min:

    :param max: Maximum
    :type max:

    :param n_ticks: Number of ticks
    :type n_ticks: int

    :param descriptor_labels: Label of descriptors
    :type descriptor_labels: list

    :param default_descriptor_pt_args: Dictionary of descriptor point arguments
    :type default_descriptor_pt_args: dict

    :param default_descriptor_label_args: Dictionary of descriptor labels
    :type default_descriptor_label_args: dict

    :param descriptor_pt_args:
    :type descriptor_pt_args: dict

    :param include_descriptors: Include the descriptors
    :type include_descriptors: bool

    :param plot_size: Size of the plot
    :type plot_size: int

    :param aspect: 
    :type aspect:

    :param subplots_adjust_kwargs: Dictionary of keyword arguments for adjusting matplotlib subplots
    :type subplots_adjust_kwargs: dict

    .. todo:: Some missing descriptions
    """
    def __init__(self):
        defaults = dict(
                resolution_enhancement = 1,
                min = None,
                max = None,
                n_ticks = 8,
                plot_function = None,
                colorbar = True,
                colormap = plt.cm.jet,
                axis_label_decimals = 2,
                log_scale = False,
                descriptor_labels = ['X_descriptor','Y_descriptor'],
                default_descriptor_pt_args = {'marker':'o'},
                default_descriptor_label_args = {},
                descriptor_pt_args = {},
                descriptor_label_args = {},
                include_descriptors = False,
                plot_size = 4,
                aspect = None,
                subplots_adjust_kwargs = {'hspace':0.35,'wspace':0.35,
                    'bottom':0.15}
                )

        for key in defaults:
            val = defaults[key]
            if not hasattr(self,key):
                setattr(self,key,val)
            elif getattr(self,key) is None:
                setattr(self,key,val)

    def update_descriptor_args(self):
        """
        Update descriptor arguments

        .. todo:: __doc__
        """
        if getattr(self,'descriptor_dict',None):
            if self.descriptor_pt_args == {}:
                for pt in self.descriptor_dict:
                    self.descriptor_pt_args[pt] = copy(
                            self.default_descriptor_pt_args)
            if self.descriptor_label_args == {}:
                for pt in self.descriptor_dict:
                    self.descriptor_label_args[pt] = copy(
                            self.default_descriptor_label_args)

    def plot_descriptor_pts(self, mapp, idx, ax, plot_in=None):
        """
        Plot descriptor points

        :param mapp:
        :type mapp:

        :param idx:
        :type idx:

        :param ax: axes object

        :param plot_in:
        :type plot_in:


        .. todo:: __doc__
        """
        if getattr(self,'descriptor_dict',None):
            self.update_descriptor_args()
            xy,rates = zip(*mapp)
            dim = len(xy[0])
            for key in self.descriptor_dict:
                pt_kwargs = self.descriptor_pt_args.get(key,
                        self.default_descriptor_pt_args)
                lab_kwargs = self.descriptor_label_args.get(key,
                        self.default_descriptor_label_args)
                if dim == 1:  # x will be descriptor values. y will be rate/coverage/etc.
                    x,y = self.descriptor_dict[key]
                    y_sp = catmap.spline(plot_in[0], plot_in[1], k=1)
                    y = y_sp(x)
                elif dim == 2:
                    x,y = self.descriptor_dict[key]
                if None not in [x,y]:
                    if pt_kwargs is not None:
                        ax.errorbar(x,y,**pt_kwargs)
                    if lab_kwargs is not None:
                        ax.annotate(key,[x,y],**lab_kwargs)
            if dim == 1:
                ax.set_xlim(self.descriptor_ranges[0])
            elif dim == 2:
                ax.set_xlim(self.descriptor_ranges[0])
                ax.set_ylim(self.descriptor_ranges[1])

    def plot_single(self, mapp, rxn_index, ax=None,
            overlay_map = None, alpha_range=None,
            **plot_args):
        """
        :param mapp:

        :param rxn_index: Index for the reaction
        :type rxn_index: int

        :param ax: axes object

        :param overlay_map:
        :type overlay_map: 

        :type alpha_range:
        :type alpha_range:

        .. todo:: __doc__
        """
        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        xy,rates = zip(*mapp)
        dim = len(xy[0])
        if dim == 1:
            x = zip(*xy)[0]
            descriptor_ranges = [[min(x),max(x)]]
            if not self.plot_function:
                if self.log_scale == True:
                    self.plot_function = 'semilogy'
                else:
                    self.plot_function = 'plot'
        elif dim == 2:
            x,y = zip(*xy)
            descriptor_ranges = [[min(x),max(x)],[min(y),max(y)]]
            if not self.plot_function:
                self.plot_function = 'contourf'
            if 'cmap' not in plot_args:
                plot_args['cmap'] = self.colormap

        eff_res =self.resolution*self.resolution_enhancement
        if self.min:
            minval = self.min
        else:
            minval = None
        maparray = RM.map_to_array(mapp,descriptor_ranges,eff_res,
                log_interpolate=self.log_scale,minval=minval)
        if self.max is None:
            self.max = maparray.T[rxn_index].max()
        if self.min is None:
            self.min = maparray.T[rxn_index].min()

        if dim == 2:
            if maparray.min() <= self.min:
                plot_args['extend'] = 'min'
            if maparray.max() >= self.max:
                plot_args['extend'] = 'max'
            if maparray.max() >= self.max and maparray.min() <= self.min:
                plot_args['extend'] = 'both'
            if 'extend' not in plot_args:
                plot_args['extend'] = 'neither'

        if self.log_scale and dim == 2:
            maparray = np.log10(maparray)
            min_val = np.log10(float(self.min))
            max_val = np.log10(float(self.max))
            if min_val < -200:
                min_val = max(maparray.min(),-200)
            elif max_val == np.inf:
                max_val = min(maparray.max(),200)
        else:
            min_val = self.min
            max_val = self.max

        maparray = np.clip(maparray,min_val,max_val)

        log_scale = self.log_scale
        if overlay_map:
            overlay_array = RM.map_to_array(overlay_map,
                    descriptor_ranges,eff_res)
            if alpha_range:
                alpha_min,alpha_max = alpha_range
            else:
                alpha_min = overlay_array.min()
                alpha_max = overlay_array.max()
            overlay_array = (overlay_array - overlay_array.min())
            overlay_array = overlay_array/(alpha_max - alpha_min)
            overlay_array = np.clip(overlay_array,0,1)
            maparray = np.clip(maparray,min_val,max_val)
            norm_array = (maparray - maparray.min())
            norm_array = norm_array/(maparray.max()-maparray.min())
            maparray = norm_array*overlay_array
            maparray = (maparray - maparray.min())
            maparray = maparray/(maparray.max()-maparray.min())
            maparray = maparray*(max_val-min_val) + min_val
            maparray=norm_array*overlay_array
            norm_array = (maparray - maparray.min())
            norm_array = norm_array/(maparray.max()-maparray.min())
            maparray = norm_array*(max_val-min_val)+min_val

        if dim == 1:
            x_range = descriptor_ranges[0]
            plot_in = [np.linspace(*x_range+eff_res),maparray[:,rxn_index]]
            plot = getattr(ax,self.plot_function)(*plot_in)
        elif dim == 2:
            x_range,y_range = descriptor_ranges
            z = maparray[:,:,rxn_index]
            if self.log_scale:
                levels = range(int(min_val),int(max_val)+1)
                if len(levels) < 3*self.n_ticks:
                    levels = np.linspace(
                            int(min_val),int(max_val),3*self.n_ticks)
            else:
                levels = np.linspace(min_val,max_val,min(eff_res,25))

            plot_in = [np.linspace(*x_range+[eff_res[0]]),
                    np.linspace(*y_range+[eff_res[1]]),z,levels]
            plot = getattr(ax,self.plot_function)(*plot_in,**plot_args)

        pos = ax.get_position()
        if self.aspect:
            ax.set_aspect(self.aspect)
            ax.apply_aspect()

        if dim == 1:
            ax.set_xlim(descriptor_ranges[0])
            ax.set_xlabel(self.descriptor_labels[0])
            ax.set_ylim([float(self.min), float(self.max)])
        elif dim == 2:
            if self.colorbar:
                if log_scale: #take only integer tick labels
                    cbar_nums = range(int(min_val),int(max_val)+1)
                    mod = max(int(len(cbar_nums)/self.n_ticks), 1)
                    cbar_nums = [n for i,n in enumerate(cbar_nums) if not i%mod]
                    cbar_nums = np.array(cbar_nums)
                else:
                    cbar_nums = np.linspace(min_val,max_val,self.n_ticks)
                formatstring = '%.'+str(self.axis_label_decimals)+'g'
                cbar_labels = [formatstring % (s,) for s in cbar_nums]
                cbar_labels = [lab.replace('e-0','e-').replace('e+0','e')
                        for lab in cbar_labels]
                plot.set_clim(min_val,max_val)
                fig = ax.get_figure()
                axpos = list(ax.get_position().bounds)
                xsize = axpos[2]*0.04
                ysize = axpos[3]
                xp = axpos[0]+axpos[2]+0.04*axpos[2]
                yp = axpos[1]
                cbar_box = [xp,yp,xsize,ysize]
                cbar_ax = fig.add_axes(cbar_box)
                cbar = fig.colorbar(mappable=plot,ticks=cbar_nums,
                        cax=cbar_ax,extend=plot_args['extend'])
                cbar.ax.set_yticklabels(cbar_labels)
                if getattr(self,'colorbar_label',None):
                    cbar_kwargs = getattr(self,'colorbar_label_kwargs',{'rotation':-90})
                    cbar_ax.set_ylabel(self.colorbar_label,**cbar_kwargs)
            if self.descriptor_labels:
                ax.set_xlabel(self.descriptor_labels[0])
                ax.set_ylabel(self.descriptor_labels[1])
            ax.set_xlim(descriptor_ranges[0])
            ax.set_ylim(descriptor_ranges[1])

        if 'title' in plot_args and plot_args['title']:
            if 'title_size' not in plot_args:
                n_pts = self.plot_size*72
                font_size = min([n_pts/len(plot_args['title']),14])
            else:
                font_size = plot_args['title_size']
            ax.set_title(plot_args['title'],size=font_size)

        if getattr(self,'n_xticks',None):
            ax.xaxis.set_major_locator(MaxNLocator(self.n_xticks))
        if getattr(self,'n_yticks',None):
            ax.yaxis.set_major_locator(MaxNLocator(self.n_yticks))

        self.plot_descriptor_pts(mapp,rxn_index,ax=ax,plot_in=plot_in)
        return ax

    def plot_separate(self,mapp,ax_list=None,indices=None,
            overlay_map = None,**plot_single_kwargs):
        """
        Generate separate plots

        .. todo:: __doc__
        """

        pts,rates = zip(*mapp)
        if indices is None:
            indices = range(0,len(rates[0]))
        n_plots = len(indices)

        if not ax_list:
            x = int(np.sqrt(n_plots))
            if x*x < n_plots:
                y = x+1
            else:
                y = x
            if x*y < n_plots:
                x = x+1
            if self.colorbar:
                fig = plt.figure(
                        figsize=(y*self.plot_size*1.25,x*self.plot_size))
            else:
                fig = plt.figure(figsize=(y*self.plot_size,x*self.plot_size))

            ax_list = []
            for i in range(0,n_plots):
                ax_list.append(fig.add_subplot(x,y,i+1))
        else:
            fig = ax_list[0].get_figure()

        if fig:
            fig.subplots_adjust(**self.subplots_adjust_kwargs)
        else:
            fig = plt.gcf()
            fig.subplots_adjust(**self.subplots_adjust_kwargs)

        plotnum = 0

        old_dict = copy(self.__dict__)

        if not self.min or not self.max:
            for id,i in enumerate(indices):
                pts, datas = zip(*mapp)
                dat_min = 1e99
                dat_max = -1e99
                for col in zip(*datas):
                    if min(col) < dat_min:
                        dat_min = min(col)
                    if max(col) > dat_max:
                        dat_max = max(col)
                if self.min is None:
                    self.min = dat_min
                if self.max is None:
                    self.max = dat_max


        for id,i in enumerate(indices):
            kwargs = plot_single_kwargs
            if self.map_plot_labels:
                try:
                    kwargs['title'] = self.map_plot_labels[i]
                except IndexError:
                    kwargs['title'] = ''

            kwargs['overlay_map'] = overlay_map
            self.__dict__.update(old_dict)
            self.plot_single(mapp,i,ax=ax_list[plotnum],**kwargs)
            plotnum+=1

        return fig

    def plot_weighted(self,mapp,ax=None,weighting='linear',
            second_map=None,indices=None,**plot_args):
        """
        Generate weighted plot

        :param mapp:
        :type mapp:

        :param ax: axes object

        :param weighting: weighting function, 'linear' or 'dual'.
        :type weighting: str

        :param second_map:

        :param indices:

        .. todo:: __doc__
        """
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            fig = ax.get_figure()

        if self.color_list is None:
            color_list = get_colors(len(mapp[0][-1])+1)
            color_list.pop(0) #remove black
        else:
            color_list = self.color_list


        pts,datas = zip(*mapp)
        if indices is None:
            indices = range(0,len(datas[0]))
        rgbs = []
        datas = zip(*datas)
        datas = [d for id,d in enumerate(datas) if id in indices]
        datas = zip(*datas)
        if second_map:
            pts2,datas2 = zip(*second_map)
            datas2 = zip(*datas2)
            datas2 = [d for id,d in enumerate(datas2) if id in indices]
            datas2 = zip(*datas2)
        else:
            datas2 = datas
        for data,data2 in zip(datas,datas2):
            if weighting=='linear':
                rs,gs,bs = zip(*color_list)
                r = 1 - sum(float((1-ri)*di) for ri,di in zip(rs,data))
                g = 1 - sum(float((1-gi)*di) for gi,di in zip(gs,data))
                b = 1 - sum(float((1-bi)*di) for bi,di in zip(bs,data))
                eff_res = self.resolution*self.resolution_enhancement
                rgbs.append([r,g,b])
            elif weighting =='dual':
                rs,gs,bs = zip(*color_list)
                r = 1 - sum(float((1-ri)*di*d2i)
                        for ri,di,d2i in zip(rs,data,data2))
                g = 1 - sum(float((1-gi)*di*d2i)
                        for gi,di,d2i in zip(gs,data,data2))
                b = 1 - sum(float((1-bi)*di*d2i)
                        for bi,di,d2i in zip(bs,data,data2))
                eff_res = 300
                rgbs.append([r,g,b])
        r,g,b = zip(*rgbs)
        x,y = zip(*pts)
        xi = np.linspace(min(x),max(x),eff_res)
        yi = np.linspace(min(y),max(y),eff_res)
        ri = griddata(x,y,r,xi,yi)
        gi = griddata(x,y,g,xi,yi)
        bi = griddata(x,y,b,xi,yi)
        rgb_array = np.zeros((eff_res,eff_res,3))
        for i in range(0,eff_res):
            for j in range(0,eff_res):
                rgb_array[i,j,0] = ri[i,j]
                rgb_array[i,j,1] = gi[i,j]
                rgb_array[i,j,2] = bi[i,j]
        xminmax,yminmax = self.descriptor_ranges
        xmin,xmax = xminmax
        ymin,ymax = yminmax
        ax.imshow(rgb_array,extent=[xmin,xmax,ymin,ymax],origin='lower')
        self.plot_descriptor_pts(mapp, i, ax)
        if getattr(self,'n_xticks',None):
            ax.xaxis.set_major_locator(MaxNLocator(self.n_xticks))
        if getattr(self,'n_yticks',None):
            ax.yaxis.set_major_locator(MaxNLocator(self.n_yticks))
        ax.set_xlabel(self.descriptor_labels[0])
        ax.set_ylabel(self.descriptor_labels[1])
        if self.aspect:
            ax.set_aspect(self.aspect)
            ax.apply_aspect()
        return fig

    def save(self, fig, save=True, default_name='map_plot.pdf'):
        """
        :param fig: figure object

        :param save: save the figure
        :type save: bool

        :param default_name: default name for the saved figure.
        :type default: str
        """
        if save == True:
            if not hasattr(self,'output_file'):
                save = default_name
            else:
                save = self.output_file
        if save:
            fig.savefig(save)

class MechanismPlot:
    """
    Class for generating potential energy diagrams

    :param energies: list of energies
    :type energies: list

    :param barriers: list of barriers
    :type barriers: list

    :param labels: list of labels
    :type labels: list
    """
    def __init__(self,energies,barriers=[],labels=[]):
        self.energies = energies
        self.barriers = barriers
        self.labels = labels
        self.energy_line_args = {'color':'k','lw':2}
        self.barrier_line_args = {'color':'k','lw':2}
        self.label_args = {'color':'k','size':16,'rotation':45}
        self.label_positions= None
        self.initial_energy = 0
        self.initial_stepnumber = 0
        self.energy_mode ='relative' #absolute
        self.energy_line_widths = 0.5

    def draw(self, ax=None):
        """
        Draw the potential energy diagram

        .. todo:: __doc__
        """
        def attr_to_list(attrname,required_length=len(self.energies)):
            """
            Return list of attributes
            :param attrname: Name of attributes
            :type attrname: list

            :param required_length: Required length for the list of attributes
            :type required_length: int

            .. todo:: __doc__
            """
            try:
                getattr(self,attrname)[0] #Ensure that it is a list
                iter(getattr(self,attrname)) #Ensure that it is a list...
                if len(getattr(self,attrname)) == required_length:
                    pass
                else:
                    raise ValueError(attrname + ' list is of length '+ \
                            str(len(getattr(self,attrname)))+ \
                            ', but needs to be of length ' + \
                            str(required_length))
                return getattr(self,attrname)
            except:
                return [getattr(self,attrname)]*required_length

        barrier_line_args = attr_to_list('barrier_line_args',
                len(self.energies)-1)
        energy_line_widths = attr_to_list('energy_line_widths')
        energy_line_args = attr_to_list('energy_line_args')
        label_args =attr_to_list('label_args')
        label_positions=attr_to_list('label_positions')

        #plot energy lines
        energy_list = np.array(self.energies)
        energy_list = (energy_list - energy_list[0])
        energy_list = list(energy_list)
        if self.energy_mode == 'relative':
            cum_energy = [energy_list[0]]
            for i,e in enumerate(energy_list[1:]):
                last = cum_energy[i]+e
                cum_energy.append(last)
            energy_list = cum_energy
        energy_list = np.array(energy_list) + self.initial_energy
        energy_list = list(energy_list)
        energy_lines = [
                [[i+self.initial_stepnumber,i+width+self.initial_stepnumber],
                    [energy_list[i]]*2]
                for i,width in enumerate(energy_line_widths)]
        self.energy_lines = energy_lines
        for i,line in enumerate(energy_lines):
            ax.plot(*line,**energy_line_args[i])

        #create barrier lines
        barrier_lines = []
        if not self.barriers: self.barriers = [0]*(len(self.energies)-1)
        for i,barrier in enumerate(self.barriers):
            xi = energy_lines[i][0][1]
            xf = energy_lines[i+1][0][0]
            yi = energy_lines[i][1][0]
            yf = energy_lines[i+1][1][0]
            if self.energy_mode == 'relative' and (barrier == 0 or barrier <= yf-yi):
                line = [[xi,xf],[yi,yf]]
                xts = (xi+xf)/2.
                yts = max([yi,yf])
            elif self.energy_mode == 'absolute' and (barrier <= yf or barrier <= yi):
                line = [[xi,xf],[yi,yf]]
                xts = (xi+xf)/2.
                yts = max([yi,yf])
            else:
                if self.energy_mode == 'relative':
                    yts = yi+barrier
                elif self.energy_mode == 'absolute':
                    yts = barrier
                    barrier = yts - yi
                barrier_rev = barrier + (yi-yf)
                if barrier > 0 and barrier_rev > 0:
                    ratio = np.sqrt(barrier)/(np.sqrt(barrier)+np.sqrt(barrier_rev))
                else:
                    print 'Warning: Encountered barrier less than 0'
                    ratio = 0.0001
                    yts = max(yi,yf)
                xts = xi + ratio*(xf-xi)
                xs = [xi,xts,xf]
                ys = [yi,yts,yf]
                f = spline(xs,ys,k=2)
                newxs = np.linspace(xi,xf,20)
                newys = f(newxs)
                line = [newxs,newys]
            barrier_lines.append(line)
        self.barrier_lines = barrier_lines
        #plot barrier lines
        for i,line in enumerate(barrier_lines):
            ax.plot(*line,**barrier_line_args[i])

        #add labels
        trans = ax.get_xaxis_transform()
        for i,label in enumerate(self.labels):
            xpos = sum(energy_lines[i][0])/len(energy_lines[i][0])
            label_position = label_positions[i]
            args = label_args[i]
            if label_position in ['top','ymax']:
                if 'ha' not in args:
                    args['ha'] = 'left'
                if 'va' not in args:
                    args['va'] = 'bottom'
                ypos = 1
                args['transform'] = trans
                ax.text(xpos,ypos,label,**args)
            elif label_position in ['bot','bottom','ymin']:
                ypos = -0.1
                ax.xaxis.set_ticks([float(sum(line[0])/len(line[0]))
                    for line in energy_lines])
                ax.set_xticklabels(self.labels)
                for attr in args.keys():
                    try:
                        [getattr(t,'set_'+attr)(args[attr])
                                for t in ax.xaxis.get_ticklabels()]
                    except:
                        pass
            elif label_position in ['omit']:
                pass
            else:
                ypos = energy_lines[i][1][0]
                if 'ha' not in args:# and 'textcoords' not in args:
                    args['ha'] = 'left'
                if 'va' not in args:# and 'textcoords' not in args:
                    args['va'] = 'bottom'
                ax.annotate(label,[xpos,ypos],**args)

class ScalingPlot:
    """
    :param descriptor_names: list of descriptor names
    :type descriptor_names: list

    :param descriptor_dict: dictionary of descriptors
    :type descriptor_dict: dict
    
    :param surface_names: list of the surface names
    :type surface_names: list 
    
    :param parameter_dict: dictionary of parameters
    :type parameter_dict: dict

    :param scaling_function: function to project descriptors into energies.
                             Should take descriptors as an argument and return a 
                             dictionary of {adsorbate:energy} pairs.
    :type scaling_function: function
    
    :param x_axis_function: function to project descriptors onto the x-axis.
                            Should take descriptors as an argument and return a
                            dictionary of {adsorbate:x_value} pairs.
    :type x_axis_function: function

    :param scaling_function_kwargs: keyword arguments for scaling_function.
    :type scaling_function_kwargs: dict

    :param x_axis_function_kwargs: keyword arguments for x_axis_function.
    :type x_axis_function_kwargs: dict
    """
    def __init__(self,descriptor_names,descriptor_dict,surface_names,
            parameter_dict,scaling_function,x_axis_function,
            scaling_function_kwargs={},x_axis_function_kwargs={},
            ):
        self.descriptor_names = descriptor_names
        self.surface_names = surface_names
        self.descriptor_dict = descriptor_dict
        self.parameter_dict = parameter_dict
        self.scaling_function = scaling_function

        self.scaling_function_kwargs = scaling_function_kwargs
        self.x_axis_function = x_axis_function

        self.x_axis_function_kwargs = x_axis_function_kwargs
        self.axis_label_size = 16
        self.surface_label_size = 16
        self.title_size = 18
        self.same_scale = True
        self.show_titles = True
        self.show_surface_labels = True
        self.subplots_adjust_kwargs = {'wspace':0.4,'hspace':0.4}
        self.x_label_dict = {}
        self.y_label_dict = {}
        self.surface_colors = []
        self.scaling_line_args = {}
        self.label_args = {}
        self.line_args = {}
        self.include_empty = True
        self.include_error_histogram = True

    def plot(self, ax_list=None, plot_size=4.0, save=None):
        """
        :param ax_list: list of axes objects
        :type ax_list: [ax]

        :param plot_size: size of the plot
        :type plot_size: float

        :param save: whether or not to save the plot
        :type save: bool

        .. todo:: __doc__
        """
        all_ads = self.adsorbate_names + self.transition_state_names
        all_ads = [a for a in all_ads if a in self.parameter_dict.keys() and
            a not in self.echem_transition_state_names]
        if self.include_empty:
            ads_names = all_ads
        else:
            ads_names = [n for n in all_ads if
                            (None in self.parameter_dict[n] or
                                sum(self.parameter_dict[n])>0.0)]

        if not self.surface_colors:
            self.surface_colors = get_colors(len(self.surface_names))
        if not self.scaling_line_args:
            self.scaling_line_args = [{'color':'k'}]*len(ads_names)
        elif hasattr(self.scaling_line_args,'update'): #its a dictionary if so.
            self.scaling_line_args = [self.scaling_line_args]*len(
                    self.adsorbate_names)
        for d in self.descriptor_names:
            if not self.include_descriptors:
                if d in ads_names:
                    ads_names.remove(d)
        if self.include_error_histogram:
            extra = 1
        else:
            extra = 0
        if not ax_list:
            spx = round(np.sqrt(len(ads_names)+extra))
            spy = round(np.sqrt(len(ads_names)+extra))
            if spy*spx < len(ads_names)+extra:
                spy+= 1
            fig = plt.figure(figsize=(spy*plot_size,spx*plot_size))
            ax_list = [fig.add_subplot(spx,spy,i+1)
                    for i in range(len(ads_names))]
        else:
            fig = None
        all_xs, all_ys = zip(*[self.descriptor_dict[s]
            for s in self.surface_names])

        fig.subplots_adjust(**self.subplots_adjust_kwargs)
        all_ys = []
        maxyrange = 0
        ymins = []
        all_err = []
        for i,ads in enumerate(ads_names):
            actual_y_vals = self.parameter_dict[ads]
            desc_vals = [self.descriptor_dict[s] for s in self.surface_names]
            scaled_x_vals = [self.x_axis_function(
                d,**self.x_axis_function_kwargs)[0][ads] for d in desc_vals]
            label = self.x_axis_function(
                desc_vals[0],**self.x_axis_function_kwargs)[-1][ads]
            scaled_y_vals = [self.scaling_function(
                d,**self.scaling_function_kwargs)[ads] for d in desc_vals]
            diffs = [scaled-actual for scaled,actual
                    in zip(scaled_y_vals,actual_y_vals) if actual != None]
            ax = ax_list[i]
            m,b = plt.polyfit(scaled_x_vals,scaled_y_vals,1)
            x_vals = np.array([round(min(scaled_x_vals),1)-0.1,
                round(max(scaled_x_vals),1)+0.1])
            ax.plot(x_vals,m*x_vals+b,**self.scaling_line_args[i])
            err = [yi - (m*xi+b) for xi,yi in zip(scaled_x_vals,actual_y_vals) if yi != None]
            all_err += err
            ax.set_xlabel(label)
            ax.set_ylabel('$E_{'+ads+'}$ [eV]')
            num_y_vals = []
#            for s,c in zip(self.surface_names,self.surface_colors):
#                print s, c
            for sf,col,x,y in zip(self.surface_names,
                    self.surface_colors,scaled_x_vals,actual_y_vals):
                if y and y != None:
                    ax.plot(x,y,'o',color=col,markersize=10,mec=None)
                    if self.show_surface_labels:
                        ax.annotate(sf,[x,y],color=col,**self.label_args)
                    num_y_vals.append(y)
            if self.show_titles:
                ax.set_title('$'+ads+'$',size=self.title_size)
            all_ys += num_y_vals
            if not num_y_vals: num_y_vals = scaled_y_vals
            dy =  max(num_y_vals) - min(num_y_vals)
            ymins.append([min(num_y_vals),max(num_y_vals)])
            if dy > maxyrange:
                maxyrange = dy
            ax.set_xlim(x_vals)
            y_range = [round(min(num_y_vals),1)-0.1,
                    round(max(num_y_vals),1)+0.1]

        self.scaling_error = all_err

        if self.same_scale == True:
            for i,ax in enumerate(ax_list):
                pad = maxyrange - (ymins[i][1]-ymins[i][0])
                y_range = [round(ymins[i][0]-pad,1)-0.1,
                        round(ymins[i][1]+pad,1)+0.1]
                ax.set_ylim(y_range)
        if self.include_error_histogram:
            err_ax = fig.add_subplot(spx,spy,len(ads_names)+1)
            err_ax.hist(all_err,bins=15)
            err_ax.set_xlabel('$E_{actual} - E_{scaled}$ [eV]')
            err_ax.set_ylabel('Counts')
            ax_list.append(err_ax)

        for ax in ax_list:
            if getattr(self,'n_xticks',None):
                ax.xaxis.set_major_locator(MaxNLocator(self.n_xticks))
            if getattr(self,'n_yticks',None):
                ax.yaxis.set_major_locator(MaxNLocator(self.n_yticks))

        if save is None:
            save = self.model_name+'_scaling.pdf'
        if save:
            fig.savefig(save)
        return fig
