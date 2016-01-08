Accessing and reformatting output
=================================

The purpose of this tutorial is to provide an explanation of how to
access the raw output of CatMAP. The plotting functions included in
CatMAP are not meant to provide publication-quality figures by default,
but rather to provide a quick way of analyzing the output and
determining whether or not the solution is correct. If you want to
create more complex plots, or more beautiful plots, then you will likely
want to reformat the raw data and plot it with software of your choice
(MATLAB, etc.) or at least know how to access the matplotlib figure
object if you are brave enough to directly edit the figure in python.

The simplest approach to post-processing is to convert the CatMAP output
into a data table and use methods you are familiar with. We will use the
CO oxidation example from tutorial 3 to provide some concrete context.
If you have run :doc:`Tutorial 3 <../tutorials/refining_a_microkinetic_model>` you should have the file "CO\_oxidation.log"
in the tutorial 3 directory. Assuming you finished the tutorial, this
file should have outputs for the coverage, rate, production\_rate, and
rate\_control. We will look at how to read in each of these and convert
them to more conventional tables.

First, lets take a look at how output is stored natively by CatMAP.
Create a script with the following commands in the tutorial 3 directory:

.. code:: python

    from catmap.model import ReactionModel

    model = ReactionModel(setup_file='CO_oxidation.log')

which reads the model (along with outputs) into the script. The model,
after its run, will have "map" attributes for each output variable. The
map is a python list structured like:

.. code:: python

    [
    [[x_0,y_0],[out_0_0, out_0_1, ... , out_0_n]],
    [[x_1,y_1],[out_1_0, out_1_1, ... , out_1_n]],
    ...,
    [[x_m,y_m],[out_m_0, out_m_1, ... , out_m_n]]
    ]

where x and y represent descriptor values, and out\_i\_j represents the
output vector at each point. For example, the coverage map on the CO
oxidation model. Add the following to the script:

.. code:: python

    for pt, cvgs in model.coverage_map:
        print 'descriptors:', pt
        print 'coverages', cvgs

If you run this you will get a bunch of numbers like:

.. code:: python

    descriptors [1.9289202008928572, 2.232142857142857]
    coverages [mpf('2.2800190488939763e-5'), mpf('0.020278026143493406'), mpf('2.2031908471388691e-7'), mpf('1.4713059814195941e-5'), mpf('0.18181715627634512')]
    descriptors [0.1428571428571428, 2.553571428571429]
    coverages [mpf('2.7893704760253176e-27'), mpf('0.050000000000000003'), mpf('2.8850535476689032e-25'), mpf('1.4466885932549904e-12'), mpf('0.94999999999855326')]
    descriptors [3.0, 1.267857142857143]
    coverages [mpf('0.049999987577682677'), mpf('4.4864510401468476e-25'), mpf('0.69752183650972651'), mpf('4.653706806446857e-10'), mpf('1.789519988373698e-17')]
    descriptors [2.1830357142857144, 2.8750000000000004]
    coverages [mpf('6.4052886383458481e-12'), mpf('0.024816133631270126'), mpf('3.4777261296001614e-13'), mpf('6.8031802914834851e-9'), mpf('0.32310023705530582')]
    ...

You should note that the coverages are "mpf" objects, which indicates
that they are multiple precision. You probably also don't know which
intermediate each coverage corresponds to. Try the following:

.. code:: python

    labels = model.output_labels['coverage']
    for pt,cvg in model.coverage_map:
        print 'descriptors',pt
        print 'intermediates',labels
        print 'coverages', [float(c) for c in cvg]

which gives the slightly more readable output:

.. code:: python

    descriptors [1.9289202008928572, 2.232142857142857]
    intermediates ('CO_s', 'O_s', 'CO_t', 'O2_t', 'O_t')
    coverages [2.2800190488939762e-05, 0.020278026143493406, 2.203190847138869e-07, 1.471305981419594e-05, 0.1818171562763451]
    descriptors [0.1428571428571428, 2.553571428571429]
    intermediates ('CO_s', 'O_s', 'CO_t', 'O2_t', 'O_t')
    coverages [2.7893704760253176e-27, 0.049999999999999996, 2.885053547668903e-25, 1.4466885932549903e-12, 0.9499999999985532]
    descriptors [3.0, 1.267857142857143]
    intermediates ('CO_s', 'O_s', 'CO_t', 'O2_t', 'O_t')
    coverages [0.049999987577682675, 4.486451040146847e-25, 0.6975218365097264, 4.653706806446857e-10, 1.7895199883736977e-17]
    descriptors [2.1830357142857144, 2.8750000000000004]
    intermediates ('CO_s', 'O_s', 'CO_t', 'O2_t', 'O_t')
    coverages [6.405288638345848e-12, 0.024816133631270124, 3.477726129600161e-13, 6.803180291483484e-09, 0.3231002370553058]
    ...

Based on this, you can probably see how to create a text table
containing coverage outputs. All the other outputs follow the same basic
format; however, there are a few tricky situations when looking at other
outputs. For example, the "labels" for reaction-specific quantities
(rates, rate constants, etc.) are actually lists which need to be
flattened into strings. Even more difficult are "matrix" outputs like
rate control, where the output is a list of lists rather than a single
list. To make life easier I have created the following script which
should create a tab-separated text table from any output (.log) file
created by CatMAP. Just place this script into the output directory, and
run it with the name of the output of interest as its first argument.

.. code:: python

    from glob import glob
    import sys
    from catmap.model import ReactionModel

    output_variable = sys.argv[1]
    logfile = glob('*.log')
    if len(logfile) > 1:
        raise InputError('Ambiguous logfile. Ensure that only one file ends with .log')
    model = ReactionModel(setup_file=logfile[0])

    if output_variable == 'rate_control':
        dim = 2
    else:
        dim = 1

    labels = model.output_labels[output_variable]

    def flatten_2d(output):
        "Helper function for flattening rate_control output"
        flat = []
        for x in output:
            flat+= x
        return flat

    #flatten rate_control labels
    if output_variable == 'rate_control':
        flat_labels = []
        for i in labels[0]:
            for j in labels[1]:
                flat_labels.append('d'+i+'/d'+j)
        labels = flat_labels

    #flatten elementary-step specific labels
    if output_variable in ['rate','rate_constant','forward_rate_constant','reverse_rate_constant']:
        str_labels = []
        for label in labels:
            states = ['+'.join(s) for s in label]
            if len(states) == 2:
                new_label = '<->'.join(states)
            else:
                new_label = states[0]+'<->'+states[1]+'->'+states[2]
            str_labels.append(new_label)
        labels = str_labels

    table = '\t'.join(list(['descriptor-'+d for d in model.descriptor_names])+list(labels))+'\n'

    for pt, output in getattr(model,output_variable+'_map'):
        if dim == 2:
            output = flatten_2d(output)
        table += '\t'.join([str(float(i)) for i in pt+output])+'\n'

    f = open(output_variable+'_table.txt','w')
    f.write(table)
    f.close()

This should give you the ability to import CatMAP output into pretty
much any other analysis or plotting program. However, if you are a
matplotlib loyalist you may want to try to edit the figure objects
directly, or perhaps even exploit the plotting capabilities of CatMAP to
plot some "map" other than those created by CatMAP. For example, lets
say that for whatever reason we wanted to plot the coverage of CO\*
times the rate of CO2 formation. We can do this by creating a python
script:

.. code:: python

    from catmap.model import ReactionModel
    from catmap.analyze import VectorMap

    log_file = 'CO_oxidation.log'
    model = ReactionModel(setup_file=log_file)

    CO_cvg_CO2_rate_map = []
    CO_idx = model.output_labels['coverage'].index('CO_s')
    CO2_idx = model.output_labels['production_rate'].index('CO2_g')

    for i,pt_cvg in enumerate(model.coverage_map):
        pt_rate = model.production_rate_map[i]
        pt,cvg = pt_cvg
        pt_i,rate = pt_rate
        assert pt == pt_i #ensure that points are the same

        CO_cvg = cvg[CO_idx]
        CO2_rate = rate[CO2_idx]
        CO_cvg_CO2_rate_map.append([pt,[CO_cvg*CO2_rate]]) #multiply the two and store in new map

    model.CO_cvg_CO2_rate_map = CO_cvg_CO2_rate_map #trick the model into thinking it has this output
    model.output_labels['CO_cvg_CO2_rate'] = ['theta_CO*r_CO2']

    vm = VectorMap(model)
    vm.plot_variable = 'CO_cvg_CO2_rate' #tell the model to plot the output you just created
    vm.log_scale = True #rates should be plotted on a log-scale
    vm.min = 1e-25 #minimum rate to plot
    vm.max = 1e3 #maximum rate to plot
    vm.threshold = 1e-25 #anything below this is considered to be 0
    vm.subplots_adjust_kwargs = {'left':0.2,'right':0.8,'bottom':0.15}
    fig = vm.plot(save='CO_cvg_CO2_rate.pdf')

If you run this script you will have a CatMAP-style plot of the CO\*
coverage multiplied by the CO2 formation rate. If you want to make
post-processing modifications to the plot, then you should note that the
output of the ``VectorMap.plot`` function is actually a
``matplotlib.figure`` object. You can get the handles for each axis by
iterating through the ``figure.axes`` attribute. Sometimes it is
convenient to label each axis the first time through to know which one
you are editing. For example, add the following lines to the script:

.. code:: python

    for j,ax in enumerate(fig.axes):
        ax.annotate(str(j), [0.05,0.9], color='w', xycoords='axes fraction')

    fig.savefig('figure_with_axes_labels.pdf')

Now if you look at the plot you will see the main axis is labeled 0,
while the colorbar is 1. You can then edit the axes properties using
matplotlib. Manipulating matplotlib figures is beyond the scope of this
tutorial, but there is plenty of good documentation at
http://matplotlib.org/.
