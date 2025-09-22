"""
AfterMD: GROMACS MD Analysis Toolkit

A comprehensive toolkit for analyzing molecular dynamics simulation results
from GROMACS using MDAnalysis and GROMACS built-in tools.
"""

__version__ = "0.1.0"
__author__ = "Research Team"

# Utils - Tool functions
from .utils import BatchProcessor, PlotManager, PathManager

# Analysis modules - Core functionality
from .analysis import (
    trajectory, structure,
    RMSDCalculator, RDFCalculator, RadiusGyrationCalculator,
    DistanceCalculator, HydrogenBondAnalyzer,
    BFactorAnalyzer, ContactMapCalculator, 
    GeometryAnalyzer, AtomInfoExtractor
)

# Preprocessing modules
from .preprocessing import PBCProcessor

# Batch processing - Simple interface
from .batch_process import process_md_tasks, discover_md_tasks, check_task_status

# SLURM cluster support
from .utils.slurm_generator import generate_slurm_scripts_for_md_tasks

__all__ = [
    # Utils
    "BatchProcessor",
    "PlotManager", 
    "PathManager",
    # Analysis submodules
    "trajectory",
    "structure", 
    # Trajectory analysis
    "RMSDCalculator",
    "RDFCalculator",
    "RadiusGyrationCalculator",
    "DistanceCalculator",
    "HydrogenBondAnalyzer",
    # Structure analysis
    "BFactorAnalyzer",
    "ContactMapCalculator",
    "GeometryAnalyzer",
    "AtomInfoExtractor",
    # Preprocessing
    "PBCProcessor",
    # Batch processing
    "process_md_tasks",
    "discover_md_tasks", 
    "check_task_status",
    # SLURM cluster support
    "generate_slurm_scripts_for_md_tasks"
]