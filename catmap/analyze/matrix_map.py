from analysis_base import *
from vector_map import *

class MatrixMap(VectorMap):
    """
    .. todo:: __doc__
    """

    def get_pts_cols(self):
        """
        .. todo:: __doc__
        """
        mapp = getattr(self,self.plot_variable+'_map')
        if not mapp:
            raise AttributeError('No output found for ' + self.plot_variable)
        pts,datas = zip(*mapp)
        cols = []
        for page in zip(*datas):
            cols += list(zip(*page))
        return pts,cols

    def include_labels_to_idxs(self):
        """
        .. todo:: __doc__
        """
        if self.include_labels:
            labels_x,labels_y = self.output_labels[self.plot_variable]
            include_x, include_y = self.include_labels
            label_idxs = []
            if not include_x: include_x = labels_x
            if not include_y: include_y = labels_y
            for x in include_x:
                x_idx = labels_x.index(x)
                for id in [labels_y.index(yi) for yi in include_y]:
                    label_idxs.append(x_idx+id)
            for y in include_y:
                y_idx = labels_y.index(y)
                for id in [labels_x.index(xi) for xi in include_x]:
                    label_idxs.append(y_idx+id*len(labels_x))

            if self.include_indices:
                print('Warning: Both labels and indices were specified.'+\
                        'The union of the sets will be used')
            else:
                self.include_indices = []
            include_indices = list(set(label_idxs+self.include_indices))
            return include_indices

    def get_labels(self):
        """
        .. todo:: __doc__
        """
        labels_x,labels_y = self.output_labels[self.plot_variable]
        if self.label_template:
            template = Template(self.label_template)
            label_strs = []
            for li in labels_x:
                for lj in labels_y:
                    label_strs.append(template.substitute({'i':li,'j':lj}))
        else:
            label_strs = []
            for li in labels_x:
                for lj in labels_y:
                    label_strs.append(str(li)+'|'+str(lj))
        return label_strs
