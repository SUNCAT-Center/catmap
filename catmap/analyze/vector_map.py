from analysis_base import *
from string import Template

class VectorMap(MapPlot, ReactionModelWrapper):
    """
    .. todo:: __doc__
    """
    def __init__(self,reaction_model):
        self._rxm = copy(reaction_model) #use copy to avoid changing model
        MapPlot.__init__(self)
        self.descriptor_labels = None
        self.log_scale = False
        self.plot_variable = 'coverage'
        self.plot_mode = 'separate'
        self.include_empty = False
        self.unique_only = False
        self.include_indices = None
        self.include_labels = None
        self.threshold = 0 #count anything less than this as 0
        self.labels = None
        self.plot_precision = 5 #round to this precision when testing for uniqueness

    def get_pts_cols(self):
        """
        .. todo:: __doc__
        """
        mapp = getattr(self, self.plot_variable+'_map')
        if not mapp:
            raise AttributeError('No output found for ' + self.plot_variable)
        pts,datas = zip(*mapp)
        # stupid hack for 1d plotting: if resolution in the y-direction is 1, make pts
        # an array of 1d vectors instead of 2d vectors.  Currently, 1d plotting can
        # only be done if resolution is set to [x_res, 1]
        if hasattr(self._rxm.resolution, '__iter__') and self._rxm.resolution[1] == 1:
            pts = [[first] for first, second in pts]
        cols = zip(*datas)
        return pts,cols

    def include_labels_to_idxs(self):
        """
        .. todo:: __doc__
        """
        if self.include_labels:
            labels = self.output_labels[self.plot_variable]
            label_idxs = [labels.index(L) for L in self.include_labels]
            if self.include_indices:
                print('Warning: Both labels and indices were specified.'+\
                        'The union of the sets will be used')
            else:
                self.include_indices = []
            include_indices = list(set(label_idxs+self.include_indices))
            return include_indices

    def get_included_indices(self,pts,cols):
        """
        .. todo:: __doc__
        """

        include_indices = self.include_labels_to_idxs()

        if not include_indices:
            include_indices = range(0,len(cols))

        if not self.include_empty:
            for i,col in enumerate(cols):
                if max([abs(c) for c in col]) <= self.threshold:
                    if i in include_indices:
                        include_indices.remove(i)

        if self.unique_only:
            def filter(val):
                val = float(val)
                if self.min:
                    val = max(val,self.min)
                if self.max:
                    val = min(val,self.max)
                if self.log_scale == True:
                    val = np.log10(abs(val))
                val = round(val,self.plot_precision)
                return val
            r_cols = [tuple([filter(ci) for ci in c]) for c in cols]
            u_cols = list(set(r_cols))
            unique_ids = [r_cols.index(c) for c in u_cols]
            include_idxs = []
            for idx in include_indices:
                if idx in unique_ids:
                    include_idxs.append(idx)
            include_indices = include_idxs
        return include_indices

    def get_labels(self):
        """
        .. todo:: __doc__
        """
        labels = self.output_labels[self.plot_variable]
        if self.label_template:
            template = Template(self.label_template)
            label_strs = [template.substitute({'i':li}) for li in labels]
        else:
            labels = [str(li) for li in labels]
        return labels

    def plot(self,ax=None,ax_list=None,save=True):
        """
        .. todo:: __doc__
        """
        pts,cols = self.get_pts_cols()
        include_indices = self.get_included_indices(pts,cols)
        datas = zip(*cols)
        mapp = zip(pts,datas)
        if not self.labels:
            self.map_plot_labels = self.get_labels()
        if not self.descriptor_labels:
            self.descriptor_labels = self.descriptor_names
        if self.plot_mode == 'separate':
            fig = self.plot_separate(mapp,ax_list,indices=include_indices)
        elif self.plot_mode == 'single':
            if ax_list:
                ax = ax_list[0]
            else:
                ax = None
            fig = self.plot_weighted(mapp,ax=ax,indices=idxs)

        self.save(fig,save=save,
                default_name = self.model_name+'_'+self.plot_variable+'.pdf')

        self._ax = ax
        self._ax_list = ax_list
        self._fig = fig
        return fig
