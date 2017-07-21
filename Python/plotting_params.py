#!/usr/bin/python
import matplotlib as mpl
from cycler import cycler

mpl.rc('lines', linewidth=2, markersize=9, markeredgewidth=1.0)
mpl.rc('font', size=22, **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)
mpl.rc('legend', framealpha=None)
mpl.rc('errorbar', capsize=4)
#mpl.rc('axes', prop_cycle=cycler(color=['#3a42d4', '#d43a3a', '#60c1dc',\
#	'#60dc96', '#dc6060', '#addc60', '#dc6079', '#60dcd4']))
# diverging
mpl.rc('axes', prop_cycle=cycler(color=['#313695', '#4575b4', '#74add1',\
	'#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026']))


