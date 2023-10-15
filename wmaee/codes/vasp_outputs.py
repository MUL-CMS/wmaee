import os
from pymatgen.io.vasp.outputs import Vasprun, Oszicar, Dos, Xdatcar
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure
from typing import Optional, Dict, Union, List
from wmaee.core.data_structs import DotDict
from monty.io import zopen
import re

def parse_output(vasp_dir: Optional[str] = None,
                 parse_dos: bool = False,
                 parse_ionic_steps: bool = False,
                 parse_outcar: bool = False,
                 return_DocDict: bool = True,
                ) -> Dict[str, Union[float, Structure, int, List[Structure], Dos, Dict]]:
    """
    Parse VASP output files in the specified directory using pymatgen.

    Args:
        vasp_dir (str, optional): The path to the directory containing VASP output files.
            If not provided, the current working directory is used.
        parse_dos (bool, optional): Whether to parse Density of States (DOS) from vasprun.xml. Defaults to False.
        parse_ionic_steps (bool, optional): Whether to parse ionic relaxation steps from XDATCAR. Defaults to False.
        parse_outcar_info (bool, optional): Whether to parse magnetization from OUTCAR. Defaults to False.
        return_DocDict (bool, optional): Whether the result to return as wmaee.core.data_structs.DocDics 
            (results.final_energy) or normal dictionary (results["final_energy"]). Defaults to True.

    Returns:
        dict: A dictionary containing parsed information from VASP outputs.
    """
    output_data = {}

    # Use the current working directory if vasp_dir is not provided
    if vasp_dir is None:
        vasp_dir = os.getcwd()

    # Check if the directory exists
    if not os.path.isdir(vasp_dir):
        raise ValueError(f"Directory '{vasp_dir}' does not exist.")

    # Parse initial and final structures from POSCAR and CONTCAR
    poscar_initial_file = os.path.join(vasp_dir, "POSCAR")
    if os.path.isfile(poscar_initial_file):
        poscar_initial = Poscar.from_file(poscar_initial_file, check_for_POTCAR=False)
        output_data["initial_structure"] = poscar_initial.structure

    poscar_final_file = os.path.join(vasp_dir, "CONTCAR")
    if os.path.isfile(poscar_final_file):
        poscar_final = Poscar.from_file(poscar_final_file, check_for_POTCAR=False)
        output_data["final_structure"] = poscar_final.structure

    # Parse vasprun.xml for electronic structure information and DOS if enabled
    vasprun_file = os.path.join(vasp_dir, "vasprun.xml")
    if parse_dos:
        if os.path.isfile(vasprun_file):
            vasprun = Vasprun(vasprun_file, parse_eigen=False, parse_potcar_file=False, parse_projected_eigen=False)
            dos = vasprun.complete_dos
            output_data["total_dos"] = dos
            output_data["final_energy"] = vasprun.final_energy # first try to read final energy

    # Parse INCAR and KPOINTS settings, and final magnetization from OUTCAR if enabled
    if parse_outcar:
        outcar_file = os.path.join(vasp_dir, "OUTCAR")
        if os.path.isfile(outcar_file):
            outcar = Outcar(outcar_file)
            output_data["total_magnetization"] = outcar.total_mag
            output_data["magnetization"] = outcar.magnetization
            output_data["final_energy"] = outcar.final_energy # second try to read final energy

    # Parse Oszicar for final energy if available
    oszicar_file = os.path.join(vasp_dir, "OSZICAR")
    if oszicar_file and os.path.isfile(oszicar_file):
        oszicar = Oszicar(oszicar_file)
        output_data["final_energy"] = oszicar.final_energy # third and final try to read final energy
        output_data["ionic_steps_energies"] = oszicar.ionic_steps
        output_data["electronic_steps_energies"] = oszicar.electronic_steps

    # Parse ionic relaxation steps from XDATCAR if enabled
    if parse_ionic_steps:
        xdatcar_file = os.path.join(vasp_dir, "XDATCAR")
        if os.path.isfile(xdatcar_file):
            xdatcar = Xdatcar(xdatcar_file)
            structures = xdatcar.structures
            output_data["ionic_relaxation_steps"] = structures

    if return_DocDict:
        return DotDict(output_data)
    else:
        return output_data
    


def parse_forces(OUTCAR: bool = False, 
                 filename: str = None, 
                 last: bool = False
                ) -> Union[List[List[List[float]]], List[List[float]]]:
    """
    Parse forces from VASP OUTCAR or vasprun.xml files.

    Args:
        OUTCAR (bool): If True, parse forces from OUTCAR; if False, parse from vasprun.xml.
        filename (str, optional): The name of the file to parse forces from. If None, uses 
            default filenames ('OUTCAR' or 'vasprun.xml').
        last (bool, optional): If True, return only forces for the last step; if False, return forces for all steps.

    Returns:
        Union[List[List[List[float]], List[List[float]]]: 
        - A list of forces for all steps from the specified file if last is False.
        - Forces for the last step from the specified file if last is True.
    """
    if OUTCAR:
        # parsing from OUTCAR
        if filename == None:
            # default filename
            filename = 'OUTCAR'
        # below definition taken of regex from https://gist.github.com/gVallverdu/0e232988f32109b5dc6202cf193a49fb
        header_pattern = r"\sPOSITION\s+TOTAL-FORCE \(eV/Angst\)\n\s-+"
        row_pattern = r"\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)"
        footer_pattern = r"\s--+"
        postprocess = lambda x: float(x)
        table_pattern_text = header_pattern + r"\s*^(?P<table_body>(?:\s+" + row_pattern + r")+)\s+" + footer_pattern
        table_pattern = re.compile(table_pattern_text, re.MULTILINE | re.DOTALL)
        rp = re.compile(row_pattern)
        # below is part of the pymatgen.io.vasp.Outcar.read_table_pattern() function
        with zopen(filename, "rt") as f:
            text = f.read()
        forces = []
        for mt in table_pattern.finditer(text):
            table_body_text = mt.group("table_body")
            table_contents = []
            for line in table_body_text.split("\n"):
                ml = rp.search(line)
                # skip empty lines
                if not ml:
                    continue
                d = ml.groupdict()
                if len(d) > 0:
                    processed_line = {k: postprocess(v) for k, v in d.items()}
                else:
                    processed_line = [postprocess(v) for v in ml.groups()]
                table_contents.append(processed_line)
            forces.append(table_contents)
        if last:
            # return only last forces
            return forces[-1]
        else:
            # return forces for all steps
            return forces
    else:
        # parsing from vasprun.xml
        if filename == None:
            # default filename
            filename = 'vasprun.xml'
        vrun = Vasprun(filename, parse_eigen=False, parse_potcar_file=False, parse_projected_eigen=False)
        tr = vrun.get_trajectory()
        if last:
            # return only last forces
            return tr[-1].site_properties['forces']
        else:
            # return forces for all steps
            return [step.site_properties['forces'] for step in tr]