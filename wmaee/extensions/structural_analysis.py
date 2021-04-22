import numpy as np

def rdf(struct, rmin=1, rmax=10, nbins=90, sc=(2,2,2)):
    """
    Calculate radial distribution function (RDF) for a given structure
    :param struct: (pymatgen.Structure or ase.Atoms) the structure
    :param rmin: (float) Rmin of the radial mesh (default: 1 (Angstroem))
    :param rmax: (float) Rmax of the radial mesh (default: 10 (Angstroem))
    :param nbins: (int) number of bins (points) along the radial mesh (default: 90)
    :param sc: ((int, int, int)) supercell size for the distance analysis (default: (2, 2, 2))
    :return: (r:numpy.array, rdf:numpy.array) radial mesh, RDF
    """

    # define radial mesh for RDF
    bin_edges = np.linspace(rmin, rmax, nbins+1)
    r = [0.5*(bin_edges[i+1]+bin_edges[i]) for i in range(len(bin_edges)-1)]
    dr = [(bin_edges[i+1]+bin_edges[i]) for i in range(len(bin_edges)-1)]

    struct.make_supercell(sc)
    dists = struct.distance_matrix.reshape(-1)
    hist, _ = np.histogram(dists, bins=bin_edges)
    return np.array(r), np.array([hist[i]/(4*np.pi*r[i]**2*dr[i]) for i in range(len(hist))])