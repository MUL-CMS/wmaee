import numpy as np
import pandas as pd
from itertools import product

def rdf(struct, rmin=1, rmax=10, nbins=90):
    """
    Calculate radial distribution function (RDF) for a given structure
    :param struct: (pymatgen.Structure or ase.Atoms) the structure
    :param rmin: (float) Rmin of the radial mesh (default: 1 (Angstroem))
    :param rmax: (float) Rmax of the radial mesh (default: 10 (Angstroem))
    :param nbins: (int) number of bins (points) along the radial mesh (default: 90)
    :return: (r:numpy.array, rdf:numpy.array) radial mesh, RDF
    """

    # define radial mesh for RDF
    bin_edges = np.linspace(rmin, rmax, nbins+1)
    r = [0.5*(bin_edges[i+1]+bin_edges[i]) for i in range(len(bin_edges)-1)]
    dr = [(bin_edges[i+1]+bin_edges[i]) for i in range(len(bin_edges)-1)]

    sc = [int(np.ceil(rmax/length)) for length in struct.lattice.abc]
    struct.make_supercell(sc)
    dists = struct.distance_matrix.reshape(-1)
    hist, _ = np.histogram(dists, bins=bin_edges)
    return np.array(r), np.array([hist[i]/(4*np.pi*r[i]**2*dr[i]) for i in range(len(hist))])

def rdf_partial(struct, rmin=1, rmax=10, nbins=90):
    """
    Calculate radial distribution function (RDF) including all partial RDFs for a given structure
    :param struct: (pymatgen.Structure or ase.Atoms) the structure
    :param rmin: (float) Rmin of the radial mesh (default: 1 (Angstroem))
    :param rmax: (float) Rmax of the radial mesh (default: 10 (Angstroem))
    :param nbins: (int) number of bins (points) along the radial mesh (default: 90)
    :return: (pandas.DataFrame) [R, RDF_tot, RDF_X-X, RDF_X-Y, ...]
    """

    bin_width = (rmax-rmin)/nbins

    # get elements for partial RDF
    species = sorted(list(struct.symbol_set))
    cols = tuple(species[i]+'-'+species[j] for i in range(len(species)) for j in range(i, len(species)))

    # create empty RDF
    RDF = pd.DataFrame(columns=('r', 'tot') + cols)
    row = {col: 0 for col in cols}
    row['tot'] = 0
    for i in range(nbins):
        row['r'] = rmin+(i+0.5)*bin_width
        RDF = RDF.append(row, ignore_index=True)
        
    # get SC sizes necessary for covering the desired range
    SC = [int(np.ceil(rmax/length)) for length in struct.lattice.abc]
    shifts = product(np.arange(-SC[0], SC[0]+1), np.arange(-SC[1], SC[1]+1), np.arange(-SC[2], SC[2]+1))
    # list of translational vertors to be applied during RDF evaluation to go beyond
    # period boundaries
    shifts = [np.dot(shift, struct.lattice.matrix) for shift in shifts]

    # for s in struct.sites:
    for s in struct.sites:
        for t in struct.sites:
            # which bond
            col = '-'.join(sorted([s.species_string, t.species_string]))
            # get all SC distances
            dist = [np.linalg.norm(s.coords-t.coords+shift) for shift in shifts]
            iR = np.floor((np.array(dist)-rmin)/bin_width).astype(int)
            for i in iR:
                if i>=0 and i<nbins:
                    RDF['tot'][i] += 1
                    RDF[col][i] += 1

    # normalize per area
    RDF['tot'] = RDF['tot']/(4*np.pi*bin_width*RDF['r']**2)
    for col in cols:
        RDF[col] = RDF[col]/(4*np.pi*bin_width*RDF['r']**2)
    
    return RDF
