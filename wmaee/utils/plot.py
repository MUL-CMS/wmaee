from typing import Optional, Tuple
import matplotlib.pyplot as plt
import numpy as np

def scatter_hist(x: np.ndarray,
                 y: np.ndarray,
                 xlabel: Optional[str] = None,
                 ylabel: Optional[str] = None,
                 figsize: Tuple[float, float] = (6, 6),
                 alpha: float = 0.1,
                 alphax: float = 0.7,
                 alphay: float = 0.7,
                 marker: str = '.') -> Tuple[plt.Figure, plt.Axes]:
    """
    Create a scatter plot with histograms on the margins.

    Parameters
    ----------
    x : np.ndarray
        Data for the x-axis.
    y : np.ndarray
        Data for the y-axis.
    xlabel : Optional[str], optional
        Label for the x-axis. Default is None.
    ylabel : Optional[str], optional
        Label for the y-axis. Default is None.
    figsize : Tuple[float, float], optional
        Size of the figure (width, height). Default is (6, 6).
    alpha : float, optional
        Transparency of the scatter plot points. Default is 0.1.
    alphax : float, optional
        Transparency of the x-axis histogram bars. Default is 0.7.
    alphay : float, optional
        Transparency of the y-axis histogram bars. Default is 0.7.
    marker : str, optional
        Marker style for the scatter plot. Default is '.'.

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
        The created figure and the main scatter plot axes.
    """

    # Create a figure with constrained layout and specified size
    fig = plt.figure(constrained_layout=True, figsize=figsize)
    
    # Create a main scatter plot axes
    ax = fig.add_gridspec(top=0.75, right=0.75).subplots()

    # Create inset axes for x-axis histogram
    ax_histx = ax.inset_axes([0, 1.05, 1, 0.25], sharex=ax)
    ax_histx.tick_params(axis='x', labelbottom=False)

    # Create inset axes for y-axis histogram
    ax_histy = ax.inset_axes([1.05, 0, 0.25, 1], sharey=ax)
    ax_histy.tick_params(axis='y', labelleft=False)

    # Plot scatter plot on the main axes with specified marker
    ax.scatter(x, y, alpha=alpha, marker=marker)

    # Plot histogram for x-axis on the inset axes
    ax_histx.hist(x, bins=100, alpha=alphax)

    # Plot histogram for y-axis on the inset axes with orientation set to horizontal
    ax_histy.hist(y, bins=100, alpha=alphay, orientation='horizontal')

    # Set labels for the x and y axes
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)

    # Return the created figure and axes
    return fig, ax