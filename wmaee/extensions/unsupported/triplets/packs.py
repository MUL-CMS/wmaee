"""
This file is used to write structure packs which then can be submitted
"""
import re
import pickle as pickle
import numpy as np
from wmaee import VASPInput, Incar, Kpoints, to_pymatgen, Potcar, vasp_interactive
from wmaee.core.types import Atoms, Collection, Optional, Union, working_directory, NoReturn
from wmaee.core.common import LoggerMixin
from typing import IO
from pprint import pprint
from os.path import exists

def create_stucture_packs(incar: Incar, structures: Collection[Atoms],
                          kpoints: Optional[Kpoints] = Kpoints.gamma_automatic((1, 1, 1)),
                          potcar: Optional[Union[Potcar, str, None]] = None, prefix: Optional[Union[None, str]] = None,
                          xc_func: Optional[Union[None, str]] = 'gga', size: Optional[int] = 500,
                          cpus: Optional[int] = 2, directory: Optional[Union[str, None, working_directory]] = None,
                          show_output: Optional[bool] = True, return_stdout: Optional[bool] = False,
                          application: Optional[Union[None, str]] = None, hostname: Optional[Union[None, str]] = None,
                          partition: Optional[Union[None, str]] = None):
    num_packs = len(structures) // size
    num_remainder = len(structures) % size
    if prefix is None:
        first_structure = to_pymatgen(next(iter(structures)))
        prefix = ''.join([s.species_string for s in first_structure])
    pack_format = '%s-%i'
    call_kwargs = dict(cpus=cpus, directory=directory, show_output=show_output, return_stdout=return_stdout,
                       application=application, hostname=hostname, partition=partition)

    with working_directory(prefix):
        sturcture_iterator = iter(structures)
        for pack_id in range(num_packs + 1):
            # prepare call_kwargs
            local_size = size if pack_id != num_packs else num_remainder
            local_offset = size * pack_id
            local_incar = incar.copy()
            local_incar['INTERACTIVE'] = True
            local_incar['NSW'] = local_size
            local_ids = [i + local_offset for i in range(local_size)]
            local_structures = [to_pymatgen(next(sturcture_iterator)) for _ in range(local_size)]
            local_directory_name = pack_format % (prefix, pack_id)
            local_call_kwargs = call_kwargs.copy()
            local_call_kwargs.update(dict(directory=local_directory_name))
            # create the VASPInput object
            local_input = VASPInput(local_structures, Incar(local_incar), kpoints=kpoints, potcar=potcar,
                                    xc_func=xc_func)
            pack_data = {
                'call_kwargs': local_call_kwargs,
                'vasp_input': local_input,
                'pack_id': pack_id,
                'structure_ids': local_ids,
                'pack_name': local_directory_name
            }
            with open('%s.pack.pickle'% local_directory_name, 'wb') as pack_handle:
                pickle.dump(pack_data, pack_handle)

def _matches(reg, line):
    return len(reg.findall(line)) > 0

class InteractiveForceExtractor(LoggerMixin):

    def __init__(self, idlist):
        super(InteractiveForceExtractor, self).__init__()
        self._idlist = idlist.copy()
        self._begin_trigger = re.compile(r'FORCES\:')
        self._end_trigger = re.compile(r'\s[0-9]+\sF\=.*')
        self._number_regex = re.compile(r'[+-]*[0]*\.[0-9]+')
        self._reading = False
        self._running_id = self._idlist.pop(0)
        self._results = {}
        self._current_buffer = []

    def write(self, line: str) -> NoReturn:
        if not self._reading:
            if _matches(self._begin_trigger, line):
                self._reading = True
        else:
            self._current_buffer.append([float(crumb) for crumb in self._number_regex.findall(line)])
            if _matches(self._end_trigger, line):
                # Now we go the final one
                forces = np.array(self._current_buffer[:-1])
                free_energy, pot_eng, diff_eng = self._current_buffer[-1]
                self._results[self._running_id] = dict(forces=forces, F=free_energy, E=pot_eng, dE=diff_eng)
                if len(self._idlist) > 0:
                    self._running_id = self._idlist.pop(0)
                else:
                    self.logger.warning('No more IDs in the ID-list. Setting all of them to undefined')
                    self._running_id = 'undefined'
                self._current_buffer = []
                self._reading = False

    @property
    def results(self):
        return self._results

def run_structure_pack(pack: Union[str, IO[bytes]]):
    if isinstance(pack, str):
        pack_handle = open(pack, 'rb')
    else:
        pack_handle = pack

    pack_data = pickle.load(pack_handle)
    pack_handle.close()
    vasp_input = pack_data['vasp_input']
    call_kwargs = pack_data['call_kwargs']
    pack_name = pack_data['pack_name']
    structure_ids = pack_data['structure_ids'].copy()
    restart_file_name = '%s.finished.pickle' % pack_name
    if exists(restart_file_name):
        with open(restart_file_name, 'rb') as finished_handle:
            finished_data = pickle.load(finished_handle)
    else:
        finished_data = {}

    finished_ids = list(finished_data.keys())
    if len(finished_ids) > 0:
        # some structures were already caluclated now lets skip them
        structure_list = list(vasp_input.structure)
        for finished_id in finished_ids:
            id_index = structure_ids.index(finished_id)
            structure_ids.pop(id_index)
            structure_list.pop(id_index)
        vasp_input.structure = structure_list
        #print('Skipping %i structures' % len(finished_ids))

    def structure_finished():
        finished_ids.append(structure_ids.pop(0))
    header = '======= %s ======='%pack_name
    print(header)
    print('='*len(header))
    pprint(call_kwargs)
    print('=' * len(header))
    force_extractor = InteractiveForceExtractor(structure_ids)
    try:
        vasp_interactive(vasp_input, **call_kwargs, callbacks=[structure_finished], output=[force_extractor])
    except RuntimeError:
        pass
    finally:
        # update the current data of finished with those new calculated
        finished_data.update(force_extractor.results)
        with open(restart_file_name, 'wb') as finished_handle:
            pickle.dump(finished_data, finished_handle)