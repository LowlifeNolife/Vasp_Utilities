"""
Add a parity checker for my sanity and make sure to test this cod eout thoroughly.
// Also add in a KPOINTS checker too!
Now just add whatever you feel like. let us now just... do the kpoints convergence setup and then be happy.
Use the functions here for the same!
"""
from pathlib import Path
from mp_api.client import MPRester

from math import gcd
from functools import reduce

from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

import vasputils.constants as const

class QueryError(Exception):
    """Custom exception for query errors."""
    pass

def obtain_material_type(mp_api_key : str,mp_id :str,semiconductor_threshold :float = 1.5) -> str:
    """
    Obtains the material type (metal or semiconductor or insulator) 
    based on the band gap value in the document for a given list of material IDs
    from the Materials Project database.

    Args:
        mp_api_key: The Materials Project API key.
        mp_id: The material ID or list of material IDs.
        semiconductor_threshold: The threshold for 
        classifying a material as a semiconductor.
    """

    with MPRester(mp_api_key) as mpr:
        try:
            docs = mpr.materials.summary.search(fields=["material_id", "formula_pretty",
                                                    "band_gap"])
        except Exception as e:
            raise QueryError(f"Error querying Materials Project: {e}")

        band_gap = docs[0].band_gap

        if band_gap == 0:
            return "metal"
        elif band_gap < semiconductor_threshold:
            return "semiconductor"
        else:
            return "insulator"

def classify_materials(mp_api_key : str, mp_ids : list[str], semiconductor_threshold : float = 1.5) -> tuple[dict]:
    """
    Classifies a list of materials as metals, semiconductors, or insulators
    based on their band gap values from the Materials Project database.

    Args:
        mp_api_key: The Materials Project API key.
        mp_ids: A list of material IDs.
        semiconductor_threshold: The threshold for classifying a material as a semiconductor.

    Returns:
        A tuple containing three dictionaries: metals, semiconductors, and insulators.

    Semimetals dependent on band structure and not the actual gap so under metals.

    simulation settings for metals same as usual.
    """
    metals = {}
    semiconductors = {}
    insulators = {}

    for mp_id in mp_ids:
        if not isinstance(mp_id, str):
            raise ValueError(f"Material ID must be a string. Got {type(mp_id)} instead.")

        material_type = obtain_material_type(mp_api_key, mp_id, semiconductor_threshold)
        if material_type == "metal":
            metals[mp_id] = material_type
        elif material_type == "semiconductor":
            semiconductors[mp_id] = material_type
        else:
            insulators[mp_id] = material_type

    return metals, semiconductors, insulators

def add_magmom_to_incar(incar : dict,magmom : int = 1) -> dict:
    """
    Adds the MAGMOM tag to the INCAR settings.

    Args:
        incar: The INCAR settings as a dictionary.
        magmom: The magnetic moment value to set for each atom.

    Returns:
        The updated INCAR settings with the MAGMOM tag added.
    """
    modify_settings = const.spin_polarisation.copy()
    modify_settings["MAGMOM"] = magmom
    return {**incar, **modify_settings}

def modify_incar_for_material_type_relaxation(incar : dict, material_type : str) -> dict:
    """
    Modifies the INCAR settings based on the material type.

    Args:
        incar: The INCAR settings as a dictionary.
        material_type: The type of material ("metal", "semiconductor", or "insulator").

    Returns:
        The updated INCAR settings based on the material type.
    """
    if material_type == "metal":
        return {**incar, **const.metal_relaxation}
    elif material_type == "semiconductor":
        return {**incar, **const.semiconductor_relaxation}
    elif material_type == "insulator":
        return {**incar, **const.insulator_relaxation}
    else:
        raise ValueError(f"Unknown material type: {material_type}")

class incar:
    def __init__(self,path : str | Path | None = None,base_dict : None | dict = const.convergence_incar_template,**kwargs):
        self.settings = base_dict.copy() if base_dict is not None else {} 
        self.settings.update(kwargs) # store the settings as a dictionary. CAN specify as many arugments.
        if path is not None:
            if isinstance(path, str):
                self.path = Path(path) / "INCAR"
            else:
                self.path = path / "INCAR"
        else:
            self.path = Path.cwd() / "INCAR"

    def write_incar(self):
        with open(self.path, 'w') as f:
            for key, value in self.settings.items():
                if isinstance(value,bool):
                    value = ".TRUE." if value else ".FALSE."
                f.write(f"{key} = {value}\n")

    def incar_update_dict(self, dict_to_add : dict): # might be redundant but useful for inplace editing.
        self.settings.update(dict_to_add)

    def incar_update_direct(self,**kwargs): 
        self.settings.update(kwargs)

    def get_incar_settings(self) -> dict:
        return self.settings

def get_kpoints_mesh_type_and_structure(mp_id : str | list[str],mp_api_key :str,MK_Pack_override : str | None = True) -> tuple[Structure,dict,str]:
        
    with MPRester(mp_api_key) as mpr:
        struct = mpr.get_structure_by_material_id(mp_id)

    real_lengths = struct.lattice.abc
    reciprocal_lengths = struct.lattice.reciprocal_lattice.abc

    real_lengths = [round(x) for x in real_lengths]
    reciprocal_lengths = [round(x) for x in reciprocal_lengths]

    abc_lcm = reduce(lambda x,y : x*y // gcd(x,y),real_lengths)
    realspace_ratio = [abc_lcm // k for k in real_lengths]

    reciprocal_gcd = reduce(lambda x,y : gcd(x,y), reciprocal_lengths)
    reciprocalspace_ratio = [k // reciprocal_gcd for k in reciprocal_lengths]

    analyzer = SpacegroupAnalyzer(struct,symprec=0.1,angle_tolerance=0.5)

    lattice_data = analyzer.get_pearson_symbol()
    key = (lattice_data[0],lattice_data[1])

    if len(const.available_meshes[key]) == 1:
        mesh_type = const.available_meshes[key][0]
    elif MK_Pack_override is True : 
        mesh_type = "Gamma"
    else:
        mesh_type = "Monkhorst-Pack"
    if len(const.available_ratios[key]) == 1:
        if const.available_ratios[key][0] == "real":
            mesh = "x".join(map(str,realspace_ratio))

        elif const.available_ratios[key][0] == "reciprocal":
            mesh ="x".join(map(str,reciprocalspace_ratio))

    else :
        mesh ="x".join(map(str,reciprocalspace_ratio))


    return struct,mesh,mesh_type