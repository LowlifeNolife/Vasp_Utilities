"""
vasputils package for setting up vasp calculations in a customisable manner.
Uses pymatgen, ase and other stuff.

Build in progess.
Things to do :
(1) Make sure POTCAR files can be customised. Can add an attribute for POTCAR pseudopotential type.
(2) How to automate this further, lets just see.
(3) By week end, want some results for band structure atleast. This one is personal

"""
import subprocess
import os
from pathlib import Path
from contextlib import contextmanager,suppress

from mp_api.client import MPRester

from pymatgen.core import Structure
from pymatgen.io.vasp import Kpoints,Poscar,Potcar
from pymatgen.io.ase import AseAtomsAdaptor
from ase import Atoms

class encut_convergence():
    def __init__(self,materials_id : (str | list[str]),
                 encut_criteria : list[int],incar_tags : dict,
                 kpoints : (str | list[str]),mp_api_key : str,kpoints_type : (str | list[str])) -> None:
        
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
        self.structures = {}
        self.kpoints = {}
        self.kpoints_mesh_type = {}
        
        if not (len(self.materials_id)== len(kpoints)== len(kpoints_type)):
           raise ValueError(
           f"Number of materials ({len(self.materials_id)}), "
           f"kpoints ({len(kpoints)}), and "
           f"kpoints_type ({len(kpoints_type)}) must be equal.")
        
        for material_id, kpoints, mesh_type in zip(self.materials_id,kpoints,kpoints_type,):
          self.kpoints[material_id] = kpoints
          self.kpoints_mesh_type[material_id] = mesh_type
    
    def get_structures(self):
        with MPRester(self.mp_api_key) as mpr:
            docs = mpr.materials.summary.search(material_ids = self.materials_id,
                                                fields = ["material_id","structure"])
            self.structures = {str(doc.material_id) : doc.structure for doc in docs}
            
    
    def construct_encut_directories(self,struct: Structure, path: Path,
                                    kpoints : str,kpoints_type : str):
        base_path = Path(path)
        for i in self.encut_criteria:
            
            encut_path = base_path / str(i)
            
            (encut_path).mkdir(parents=True,exist_ok=True)
            
            file_path = encut_path / "INCAR"
            
            with open(file_path,"w") as o:   # INCAR CUSTOMISED DUE TO Incar not working properly
                incar = {**self.incar_tags , "ENCUT" : i}
                for name,type in incar.items():
                    o.write(f"{name} = {type}\n")
            
            kpoints_path =  encut_path / "KPOINTS"
            kpts = list(map(int,kpoints.split()))
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
            ptcr = Potcar(symbols=symbols,functional="PBE_54")
            ptcr.write_file(potcar_path)  
 
    def setup_material_encut_directories(self,path :Path | str |None = None):
    
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
                                            self.kpoints[id],self.kpoints_mesh_type[id])
    
        