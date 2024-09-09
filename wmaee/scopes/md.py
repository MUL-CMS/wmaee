"""Analysis and plotting of (AI)MD runs"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional, Tuple, List

def plot_MD(df: pd.DataFrame, 
            time: Optional[str] = None,
            timestep: int = 1,
            grid: Tuple[int, int] = (1, 1),
            props: List[str] = [],
            labels: Optional[List[str]] = None,
            show: bool = True,
            **kwargs) -> Tuple[plt.Figure, np.ndarray]:
    """
    Plots multiple properties from a DataFrame in a grid layout.
    
    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame containing the data to plot.
    time : Optional[str], optional
        The column name to use as the time axis. If None, a default range is used, by default None.
    timestep : int, optional
        The timestep interval for the default time axis, by default 1.
    grid : Tuple[int, int], optional
        The grid size for the subplots (rows, columns), by default (1, 1).
    props : List[str], optional
        The list of DataFrame columns to plot, by default [].
    labels : Optional[List[str]], optional
        The labels for the y-axis of each subplot, by default None.
    show : bool, optional
        Whether to display the plot, by default True.
    **kwargs
        Additional keyword arguments passed to `plt.subplots`.

    Returns
    -------
    Tuple[plt.Figure, np.ndarray]
        The figure and axes array of the created subplots.
    """
    
    # Determine the time axis values
    if time is None:
        time = np.arange(len(df)) * timestep
    else:
        time = df[time]    
    
    # Create subplots with the specified grid size
    fig, axs = plt.subplots(*grid, squeeze=False, sharex=True, **kwargs)
            
    # Iterate through the grid and plot the data
    for i in range(grid[0]):
        for j in range(grid[1]):
            n = i * grid[1] + j
            if n >= len(props):
                break
            # Plot the data for the nth property
            axs[i, j].plot(time, df[props[n]])
            # Set the y-axis label if provided
            if labels is not None:
                axs[i, j].set_ylabel(labels[n])
            # Set the x-axis label for the last row
            axs[-1, j].set_xlabel('step')
        else:
            continue
        break

    # Adjust subplot spacing
    fig.subplots_adjust(hspace=0)
    # Show the plot if requested
    if show:
        plt.show()
    
    return fig, axs