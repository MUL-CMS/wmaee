from pyiron_atomistics.atomistics.job.atomistic import AtomisticGenericJob
from pyiron_atomistics import Atoms, pyiron_to_ase
from pyiron_base import GenericParameters
from ase.io import write
from os.path import join
import pandas as pd


class CHGNet(AtomisticGenericJob):
    def __init__(self, project, job_name):
        super(CHGNet, self).__init__(project, job_name)
        self.structure = None 
        self.input = CalcParams()
        self._executable_activate(codename='chgnet')

    @property
    def fmax(self):
        return self.input.get('fmax')
            
    @fmax.setter
    def fmax(self, fmax: float):
        self.input.set(fmax=fmax)
        
    
    def write_input(self):
        write(
            join(self.working_directory, 'structure.cif'), 
            pyiron_to_ase(self.structure)
        )
        if self._generic_input['calc_mode'] == 'static':
            self._write_calc_static()

            
    def _write_calc_static(self):
        params = []
        for p in self.input.keys():
            params.append(f'{p}={self.input.get(p)}')
        params = ','.join(params)
        script = [
            'from chgnet.model import CHGNet, StructOptimizer',
            'from ase.io import read',
            'struct = read("structure.cif")',
            'relaxer = StructOptimizer()',            
            'result = relaxer.relax(struct, save_path="relax.out",'+params+')'
        ]
        with open(join(self.working_directory, 'calc_script.py'), 'w') as f:
            f.writelines("\n".join(script))

            
    def _parse_calc_static(self, output='error.out'):
        try:
            header = pd.read_csv(
                join(self.working_directory, output),
                skiprows=2, 
                nrows=0,
                sep='\s+',
            )
            df = pd.read_csv(
                join(self.working_directory, output), 
                skiprows=3, 
                sep='\s+',
                header=None
            )
            df.columns = ['optimizer_class'] + list(header.columns[:])
            return df
        except:
            self.logger.warning(f'Cannot partse output from {output} file.')
        
        
        
    def collect_output(self):
        if self._generic_input['calc_mode']=='static':
            df = self._parse_calc_static()
            with self.project_hdf5.open('output/generic') as h5out:
                h5out['energy_pot'] = df['Energy'].values
                h5out['energy_tot'] = df['Energy'].values
                h5out['max_force'] = df['fmax'].values

                

class CalcParams(GenericParameters):
    def __init__(self, input_file_name=None, table_name="calc_params"):
        super(CalcParams, self).__init__(
            table_name=table_name,
        )
    
