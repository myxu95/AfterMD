from . import trajectory
from . import structure

# Import individual modules for convenience
from .trajectory import (
    RMSDCalculator,
    RDFCalculator,
    RadiusGyrationCalculator,
    DistanceCalculator,
    HydrogenBondAnalyzer
)

from .structure import (
    BFactorAnalyzer,
    ContactMapCalculator,
    GeometryAnalyzer,
    AtomInfoExtractor
)

__all__ = [
    # Submodules
    "trajectory",
    "structure",
    # Trajectory analysis modules
    "RMSDCalculator",
    "RDFCalculator",
    "RadiusGyrationCalculator",
    "DistanceCalculator",
    "HydrogenBondAnalyzer",
    # Structure analysis modules
    "BFactorAnalyzer",
    "ContactMapCalculator",
    "GeometryAnalyzer",
    "AtomInfoExtractor"
]