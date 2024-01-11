"""
Collection of routines for reading and plotting density of states.

David Holec
david.holec@unileoben.ac.at
"""

import numpy as np
from wmaee.core.requirements import test_pmg


class dos():
    
    def __init__(self):
        """
        Creates new instance of DOS.
        """
        
        self.energy = {'Emin': None, 'Emax': None, 'EF': None, 'EN': None, 'Eg': None, 'E': []}
        self.TDOS = np.array([])
        self.PDOS: None
        self.struct: None



    # @TODO spin-polarized DOS!!!
    def read_DOS_VASP(self, fDOSCAR='DOSCAR', fPOSCAR='POSCAR', ase_atoms=True):
        """
        Reads DOSCAR file (VASP output).
        
        Arguments:
            fDOSCAR:  File with DOSCAR format.
                      Default: DOSCAR
            fPOSCAR:  File with POSCAR format (for assignment of PDOS).
                      Default: POSCAR            
        """
        
        if not ase_atoms and test_pmg():
            from pymatgen.io.vasp import Poscar
            self.struct = Poscar.from_file(fPOSCAR).structure
            species = [s.species_string for s in self.struct.sites]
        else:
            from ase.io.vasp import read_vasp
            self.struct = read_vasp(fPOSCAR)
            species = self.struct.get_chemical_symbols()
        
        with open(fDOSCAR) as f:
            DOSCAR = f.readlines()
        
        self.energy['Emin'] = float(DOSCAR[5].split()[1])
        self.energy['Emax'] = float(DOSCAR[5].split()[0])
        self.energy['EN'] = int(DOSCAR[5].split()[2])
        self.energy['EF'] = float(DOSCAR[5].split()[3])
        
        E = []
        TDOS = []
        for line in range(6, 6+self.energy['EN']):
            E.append(float(DOSCAR[line].split()[0]))
            TDOS.append(float(DOSCAR[line].split()[1]))
            
        self.energy['E'] = np.array(E)
        self.TDOS = np.array(TDOS)
        self.energy['Eg'] = 0
        for ene, dos in zip(self.energy['E'], self.TDOS):
            if ene > self.energy['EF'] and dos < 1e-5:
                self.energy['Eg'] = ene - self.energy['EF']
            elif self.energy['Eg'] > 0 and dos > 1e-5:
                break
        
        if(len(DOSCAR) == 5+(len(species)+1)*(self.energy['EN']+1)):
            # print('info: PDOS seems to be available')
            self.PDOS = []
            for i in range(len(species)):
                LDOS = []
                s, p, d = [], [], []
                px, py, pz = [], [], []
                dxy, dyz, dz2, dxz, dx2y2 = [], [], [], [], []
                for line in range(6+(i+1)*(self.energy['EN']+1), 6+(i+2)*(self.energy['EN']+1)-1):
                    pdos = [float(i) for i in DOSCAR[line].split()]
                    LDOS.append(sum(pdos[1:]))
                    s.append(pdos[1]), p.append(sum(pdos[2:5])), d.append(sum(pdos[5:10]))
                    px.append(pdos[2]), py.append(pdos[3]), pz.append(pdos[4])
                    dxy.append(pdos[5]), dyz.append(pdos[6]), dz2.append(pdos[7]), dxz.append(pdos[8]), dx2y2.append(pdos[9])        
                self.PDOS.append({'species': species[i],
                                  'LDOS': np.array(LDOS),
                                  's': np.array(s), 'p': np.array(p), 'd': np.array(d),
                                  'px': np.array(px), 'py': np.array(py), 'pz': np.array(pz),
                                  'dxy': np.array(dxy), 'dyz': np.array(dyz), 'dz2': np.array(dz2), 
                                  'dxz': np.array(dxz), 'dx2y2': np.array(dx2y2)})
        #self.DOS = DOS
