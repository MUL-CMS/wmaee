"""Analysis and plotting of (AI)MD runs"""

import matplotlib.pyplot as plt
import numpy as np
from typing import Optional

def plot_MD(df, 
            time: Optional[str] = None,
            timestep = 1,
            grid = (1, 1),
            props = [],
            show = True,
            **kwargs):
    
    if time == None:
        x = np.array(range(len(df)))*timestep
        
    fig, ax = plt.subplots(*grid, **kwargs)
    ax[0][0].plot(x, df[props[0]])
    
    fig.show()