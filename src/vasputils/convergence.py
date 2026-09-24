"""
vasputils package for setting up vasp calculations in a customisable manner.
Uses pymatgen, ase and other stuff.

Build in progess.
Things to do :
(1) Make sure POTCAR files can be customised. Can add an attribute for POTCAR pseudopotential type.
(2) How to automate this further, lets just see.
(3) By week end, want some results for band structure atleast. This one is personal
(4) Add a preprocessing step

Note : Am adding a POTCAR map as monkey wrench fix 
"""
POTCAR_MAP = {
    "Nb": "Nb_pv",
}

import subprocess
import os
from pathlib import Path
from contextlib import contextmanager,suppress

from mp_api.client import MPRester

from pymatgen.core import Structure
from pymatgen.io.vasp import Kpoints,Poscar,Potcar
from pymatgen.io.ase import AseAtomsAdaptor
from ase import Atoms
from preprocessing import incar,get_kpoints_mesh_type_and_structure,modify_incar_for_material_type_relaxation,obtain_material_type
class encut_convergence():
    def __init__(self,materials_id : (str | list[str]),
                 encut_criteria : list[int],incar_tags : dict | None,
                 mp_api_key : str,MK_Pack_override: bool = True) -> None:
        
        if isinstance(kpoints,str):
            kpoints = [kpoints]
            
        if isinstance(kpoints_type,str):
            kpoints_type = [kpoints_type]
        
        if isinstance(materials_id,str):
            materials_id = [materials_id]
            
        self.materials_id = materials_id
        self.encut_criteria = encut_criteria
        self.incar_tags = incar_tags
        self.mp_api_key = mp_api_key
        self.MK_Pack_override = MK_Pack_override
        
        self.structures = {}
        self.kpoints = {}
        self.kpoints_mesh_type = {}
        self.material_types = {}
        
        

    def get_material_data(self):
        for material_id in self.materials_id:
            struct,mesh,mesh_type = get_kpoints_mesh_type_and_structure(material_id,
                                                                        self.mp_api_key,
                                                                        self.MK_Pack_override
                                                                        )
            
            self.structures[material_id] = struct
            self.kpoints[material_id] = mesh
            self.kpoints_mesh_type[material_id] = mesh_type
        
        for material_id in self.materials_id:
            type = obtain_material_type(self.mp_api_key,
                                        mp_id=material_id)
            
            self.material_types[material_id] = type
    
            
    def construct_encut_directories(self,struct: Structure, path: Path,
                                    kpoints : str,kpoints_type : str,material_type : str):
        base_path = Path(path)
        for i in self.encut_criteria:
            
            encut_path = base_path / str(i)
            
            (encut_path).mkdir(parents=True,exist_ok=True)
            
            incar = modify_incar_for_material_type_relaxation(incar=self.incar_tags,material_type=material_type)
            
            file_path = encut_path / "INCAR"
            
            incar["ENCUT"] = i
            
            with open(file_path,"w") as o:   # INCAR CUSTOMISED DUE TO Incar not working properly
                
                for name,type in incar.items():
                    o.write(f"{name} = {type}\n")
            
            kpoints_path =  encut_path / "KPOINTS"
            kpts = list(map(int,kpoints.split("x")))
            kp = Kpoints(comment = "Kpoints for this ENCUT convergence",
                         style = kpoints_type,kpts = [kpts],
                         kpts_shift=(0,0,0))
            
            kp.write_file(kpoints_path)
            
            poscar_path = encut_path / "POSCAR"
            
            pscr = Poscar(structure=struct,
                          comment = f"POSCAR file for {struct.composition.reduced_formula}")
            pscr.write_file(poscar_path)
            
            potcar_path = encut_path / "POTCAR"
            symbols = pscr.site_symbols
            symbols = [POTCAR_MAP.get(symbol, symbol) for symbol in pscr.site_symbols]
            
            ptcr = Potcar(symbols=symbols,functional="PBE_54")
            ptcr.write_file(potcar_path)  
 
    def setup_material_encut_directories(self,path :Path | str | None = None):
    
        if path is None: 
           path = Path.cwd()
        elif type(path) == str:
            path = Path(path)
        print("Structures:")
        print(self.structures.keys())

        print("Kpoints:")
        print(self.kpoints.keys())
        
        print("Kpoint mesh types:")
        print(self.kpoints_mesh_type.keys())
        for id,struct in self.structures.items() :
           name = struct.composition.reduced_formula
           base_path = path /name
           (base_path).mkdir(parents=True,exist_ok = True)
           self.construct_encut_directories(struct,base_path,
                                            self.kpoints[id],self.kpoints_mesh_type[id],self.material_types[id])
    
class kpoints_convergence():
    def __init__():
        pass
    