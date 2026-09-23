# ---------------------------------------------------------------------------
# Fixed mappings for different uses.
# ---------------------------------------------------------------------------
# Pearson Symbol Components in Crystallography

crystal_systems = {
    "a": "Triclinic",
    "m": "Monoclinic",
    "o": "Orthorhombic",
    "t": "Tetragonal",
    "h": "Hexagonal (and Trigonal)",
    "c": "Cubic"
}

lattice_variants = {
    "P": "Primitive",
    "I": "Body-centered",
    "F": "Face-centered",
    "C": "Side-centered (A-, B-, or C-face-centered)",
    "R": "Rhombohedral"
}

available_meshes = {
        ("a", "P"): ["Gamma", "Monkhorst-Pack"],

        ("m", "P"): ["Gamma", "Monkhorst-Pack"],
        ("m", "C"): ["Gamma", "Monkhorst-Pack"],

        ("o", "P"): ["Gamma", "Monkhorst-Pack"],
        ("o", "C"): ["Gamma", "Monkhorst-Pack"],
        ("o", "I"): ["Gamma", "Monkhorst-Pack"],
        ("o", "F"): ["Gamma"],

        ("t", "P"): ["Gamma", "Monkhorst-Pack"],
        ("t", "I"): ["Gamma", "Monkhorst-Pack"],

        ("h", "R"): ["Gamma"],
        ("h", "P"): ["Gamma"],

        ("c", "P"): ["Gamma", "Monkhorst-Pack"],
        ("c", "I"): ["Gamma", "Monkhorst-Pack"],
        ("c", "F"): ["Gamma"]
    }
available_ratios = {
        ("a", "P"): ["real", "reciprocal"],

        ("m", "P"): ["real", "reciprocal"],
        ("m", "C"): ["real", "reciprocal"],

        ("o", "P"): ["real", "reciprocal"],
        ("o", "C"): ["real", "reciprocal"],
        ("o", "I"): ["real"],
        ("o", "F"): ["reciprocal"],

        ("t", "P"): ["real", "reciprocal"],
        ("t", "I"): ["real"],

        ("h", "R"): ["real", "reciprocal"],
        ("h", "P"): ["real", "reciprocal"],

        ("c", "P"): ["real", "reciprocal"],
        ("c", "I"): ["real", "reciprocal"],
        ("c", "F"): ["real", "reciprocal"]
}



from copy import deepcopy
from contextlib import contextmanager
from enum import Enum
from typing import Generator


# ---------------------------------------------------------------------------
# Default INCAR templates
# ---------------------------------------------------------------------------

metal_relaxation = {
    "ISMEAR": 1,
    "SIGMA": 0.2,
}

metal_relaxation_accurate = {
    "ISMEAR": 2,
    "SIGMA": 0.2,
}

semiconductor_relaxation = {
    "ISMEAR": 0,
    "SIGMA": 0.05,
}

insulator_relaxation = {
    "ISMEAR": 0,
    "SIGMA": 0.01,
}

density_of_states_3D = {
    "ISMEAR": -5,
}

generic_relaxation = {
    "ISMEAR": 0,
    "SIGMA": 0.05,
}

density_of_states_2D = {
    # Tetrahedron method requires a 3D mesh.
    # For 2D materials we use Gaussian smearing.
    "ISMEAR": 0,
    "SIGMA": 0.02,
}

spin_polarisation = {
    "ISPIN": 2,
    # Add MAGMOM separately according to the system.
}

convergence_incar_template = {
    # System & Precision Info
    "SYSTEM": "Unknown_Material_Convergence_Test",
    "PREC": "Accurate",
    "LASPH": True,

    # Static electronic run
    "IBRION": -1,
    "NSW": 0,

    # Electronic convergence
    "ALGO": "Normal",
    "NELM": 150,
    "EDIFF": 1e-6,

    # Smearing
    "ISMEAR": 0,
    "SIGMA": 0.05,

    # Performance & minimal disk I/O
    "LREAL": False,
    "LWAVE": False,
    "LCHARG": False,
}


# ---------------------------------------------------------------------------
# Strict list of available templates. Add more if required.
# ---------------------------------------------------------------------------

class INCARTemplate(Enum):
    METAL = "metal_relaxation"
    METAL_ACCURATE = "metal_relaxation_accurate"
    SEMICONDUCTOR = "semiconductor_relaxation"
    INSULATOR = "insulator_relaxation"
    DOS_3D = "density_of_states_3D"
    DOS_2D = "density_of_states_2D"
    GENERIC_RELAXATION = "generic_relaxation"
    CONVERGENCE = "convergence_incar_template"
    SPIN_POLARISATION = "spin_polarisation"


# ---------------------------------------------------------------------------
# Mapping between Enum members and actual INCAR dictionaries. Add the necessary mappings for other stuff too.
# ---------------------------------------------------------------------------

INCAR_TEMPLATES = {
    INCARTemplate.METAL: metal_relaxation,
    INCARTemplate.METAL_ACCURATE: metal_relaxation_accurate,
    INCARTemplate.SEMICONDUCTOR: semiconductor_relaxation,
    INCARTemplate.INSULATOR: insulator_relaxation,
    INCARTemplate.DOS_3D: density_of_states_3D,
    INCARTemplate.DOS_2D: density_of_states_2D,
    INCARTemplate.GENERIC_RELAXATION: generic_relaxation,
    INCARTemplate.CONVERGENCE: convergence_incar_template,
    INCARTemplate.SPIN_POLARISATION: spin_polarisation,
}


# ---------------------------------------------------------------------------
# Preserve the original defaults
# ---------------------------------------------------------------------------

_DEFAULTS = {
    template: deepcopy(settings)
    for template, settings in INCAR_TEMPLATES.items()
}


# ---------------------------------------------------------------------------
# Functions for working with the templates
# ---------------------------------------------------------------------------

def get(
    template: INCARTemplate,
) -> dict:
    """
    Return a copy of an INCAR template.

    The returned dictionary can be modified without changing
    the stored template.
    """

    if not isinstance(template, INCARTemplate):
        raise TypeError(
            f"template must be an INCARTemplate, "
            f"got {type(template).__name__}"
        )

    return deepcopy(INCAR_TEMPLATES[template])


def modify(
    template: INCARTemplate,
    **kwargs,
) -> None:
    """
    Permanently modify an INCAR template for the current
    Python session.

    Example:
        modify(INCARTemplate.METAL, ISMEAR=2, SIGMA=0.1)
    """

    if not isinstance(template, INCARTemplate):
        raise TypeError(
            f"template must be an INCARTemplate, "
            f"got {type(template).__name__}"
        )

    INCAR_TEMPLATES[template].update(kwargs)


def reset(
    template: INCARTemplate,
) -> None:
    """
    Restore an INCAR template to its original default values.
    """

    if not isinstance(template, INCARTemplate):
        raise TypeError(
            f"template must be an INCARTemplate, "
            f"got {type(template).__name__}"
        )

    target = INCAR_TEMPLATES[template]
    default = _DEFAULTS[template]

    target.clear()
    target.update(deepcopy(default))


@contextmanager
def temporary(
    template: INCARTemplate,
    **kwargs,
) -> Generator[dict,None,None]:
    """
    Temporarily modify an INCAR template.

    The original values are restored when the context exits,
    even if an exception occurs.

    Example:

        with temporary(INCARTemplate.METAL, ISMEAR=2):
            ...
    """

    if not isinstance(template, INCARTemplate):
        raise TypeError(
            f"template must be an INCARTemplate, "
            f"got {type(template).__name__}"
        )

    target = INCAR_TEMPLATES[template]

    old_values = deepcopy(target)

    try:
        target.update(kwargs)
        yield target

    finally:
        target.clear()
        target.update(old_values)