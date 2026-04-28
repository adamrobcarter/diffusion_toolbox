
import matplotlib

# matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=["#5BCEFA", "#F5A9B8"]) 
discrete_colors = ["#00b3fa", "#ff869e", "#736B92", 'tab:red']
matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=discrete_colors)

# https://duetosymmetry.com/code/latex-mpl-fig-tips/
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{siunitx-v2}' # sisetup detect all
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = 'Computer Modern'
matplotlib.rcParams['font.size'] = 8

# \textwidth is currently 345pt
pt = 1/72 # 72pt = 1in
textwidth = 345*pt
figsize_small = (3, 2.7)
figsize_halfpage = (textwidth/2, 2.2)