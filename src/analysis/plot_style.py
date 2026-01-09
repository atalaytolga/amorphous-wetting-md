# src/analysis/plot_style.py

import matplotlib.pyplot as plt
import matplotlib as mpl

def set_plot_style():
    """
    Applies style for thesis plots.
    """
    mpl.rcParams.update({
        # FONT
        'font.family': 'serif',
        'font.serif': ['DejaVu Serif', 'Times New Roman', 'serif'],
        'mathtext.fontset': 'cm', 
        'font.size': 12,
        
        # AXES & TICKS
        'axes.labelsize': 14,
        'axes.titlesize': 14,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.top': True,
        'ytick.right': True,
        'xtick.major.size': 6,
        'ytick.major.size': 6,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11,
        
        # LINES & MARKERS 
        'lines.linewidth': 1.0,     
        'lines.markersize': 4,      
        'lines.markeredgewidth': 0, 
        
        # FIGURE
        'figure.figsize': (6, 4.5), 
        'figure.dpi': 120,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        
        # GRID & LEGEND
        'axes.grid': True,
        'grid.alpha': 0.3,
        'grid.linestyle': '--',
        'legend.frameon': False,    
        'legend.fontsize': 11
    })

