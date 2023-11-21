# Patching pymatgen Vasprun class to read in also stresses from the MD run

from wmaee.core.requirements import test_pmg
if test_pmg():

    from pymatgen.io.vasp.outputs import *
    try:
        from pymatgen.io.vasp.outputs import _parse_varray
    except:
        # the function has been renamed in v2023.11.12
        from pymatgen.io.vasp.outputs import _parse_vasp_array as _parse_varray

    # Added lines: 105--106
    def _parse_patch(self, stream, parse_dos, parse_eigen, parse_projected_eigen):
            self.efermi = self.eigenvalues = self.projected_eigenvalues = self.projected_magnetisation = None
            self.dielectric_data = {}
            self.other_dielectric = {}
            self.incar = {}
            ionic_steps = []

            md_data = []
            parsed_header = False
            try:
                for _, elem in ET.iterparse(stream):
                    tag = elem.tag
                    if not parsed_header:
                        if tag == "generator":
                            self.generator = self._parse_params(elem)
                        elif tag == "incar":
                            self.incar = self._parse_params(elem)
                        elif tag == "kpoints":
                            if not hasattr(self, "kpoints"):
                                self.kpoints, self.actual_kpoints, self.actual_kpoints_weights = self._parse_kpoints(elem)
                        elif tag == "parameters":
                            self.parameters = self._parse_params(elem)
                        elif tag == "structure" and elem.attrib.get("name") == "initialpos":
                            self.initial_structure = self._parse_structure(elem)
                            self.final_structure = self.initial_structure
                        elif tag == "atominfo":
                            self.atomic_symbols, self.potcar_symbols = self._parse_atominfo(elem)
                            self.potcar_spec = [{"titel": p, "hash": None} for p in self.potcar_symbols]
                    if tag == "calculation":
                        parsed_header = True
                        if not self.parameters.get("LCHIMAG", False):
                            ionic_steps.append(self._parse_calculation(elem))
                        else:
                            ionic_steps.extend(self._parse_chemical_shielding_calculation(elem))
                    elif parse_dos and tag == "dos":
                        try:
                            self.tdos, self.idos, self.pdos = self._parse_dos(elem)
                            self.efermi = self.tdos.efermi
                            self.dos_has_errors = False
                        except Exception:
                            self.dos_has_errors = True
                    elif parse_eigen and tag == "eigenvalues":
                        self.eigenvalues = self._parse_eigen(elem)
                    elif parse_projected_eigen and tag == "projected":
                        self.projected_eigenvalues, self.projected_magnetisation = self._parse_projected_eigen(elem)
                    elif tag == "dielectricfunction":
                        if (
                            "comment" not in elem.attrib
                            or elem.attrib["comment"] == "INVERSE MACROSCOPIC DIELECTRIC TENSOR (including "
                            "local field effects in RPA (Hartree))"
                        ):
                            if "density" not in self.dielectric_data:
                                self.dielectric_data["density"] = self._parse_diel(elem)
                            elif "velocity" not in self.dielectric_data:
                                # "velocity-velocity" is also named
                                # "current-current" in OUTCAR
                                self.dielectric_data["velocity"] = self._parse_diel(elem)
                            else:
                                raise NotImplementedError("This vasprun.xml has >2 unlabelled dielectric functions")
                        else:
                            comment = elem.attrib["comment"]
                            # VASP 6+ has labels for the density and current
                            # derived dielectric constants
                            if comment == "density-density":
                                self.dielectric_data["density"] = self._parse_diel(elem)
                            elif comment == "current-current":
                                self.dielectric_data["velocity"] = self._parse_diel(elem)
                            else:
                                self.other_dielectric[comment] = self._parse_diel(elem)

                    elif tag == "varray" and elem.attrib.get("name") == "opticaltransitions":
                        self.optical_transition = np.array(_parse_varray(elem))
                    elif tag == "structure" and elem.attrib.get("name") == "finalpos":
                        self.final_structure = self._parse_structure(elem)
                    elif tag == "dynmat":
                        hessian, eigenvalues, eigenvectors = self._parse_dynmat(elem)
                        # n_atoms is not the total number of atoms, only those for which force constants were calculated
                        # https://github.com/materialsproject/pymatgen/issues/3084
                        n_atoms = len(hessian) // 3
                        hessian = np.array(hessian)
                        self.force_constants = np.zeros((n_atoms, n_atoms, 3, 3), dtype="double")
                        for ii in range(n_atoms):
                            for jj in range(n_atoms):
                                self.force_constants[ii, jj] = hessian[ii * 3 : (ii + 1) * 3, jj * 3 : (jj + 1) * 3]
                        phonon_eigenvectors = []
                        for ev in eigenvectors:
                            phonon_eigenvectors.append(np.array(ev).reshape(n_atoms, 3))
                        self.normalmode_eigenvals = np.array(eigenvalues)
                        self.normalmode_eigenvecs = np.array(phonon_eigenvectors)
                    elif self.incar.get("ML_LMLFF"):
                        if tag == "structure" and elem.attrib.get("name") is None:
                            md_data.append({})
                            md_data[-1]["structure"] = self._parse_structure(elem)
                        elif tag == "varray" and elem.attrib.get("name") == "forces":
                            md_data[-1]["forces"] = _parse_varray(elem)
    # +++++++++++++++++++++++++++++++++
                        elif tag == "varray" and elem.attrib.get("name") == "stress":
                            md_data[-1]["stress"] = _parse_varray(elem)
    # +++++++++++++++++++++++++++++++++
                        elif tag == "energy":
                            d = {i.attrib["name"]: float(i.text) for i in elem.findall("i")}
                            if "kinetic" in d:
                                md_data[-1]["energy"] = {i.attrib["name"]: float(i.text) for i in elem.findall("i")}
            except ET.ParseError as exc:
                if self.exception_on_bad_xml:
                    raise exc
                warnings.warn(
                    "XML is malformed. Parsing has stopped but partial data is available.",
                    UserWarning,
                )
            self.ionic_steps = ionic_steps
            self.md_data = md_data
            self.vasp_version = self.generator["version"]
            
    Vasprun._parse = _parse_patch