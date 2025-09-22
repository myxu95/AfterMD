#!/usr/bin/env python3
"""
Basic usage examples for AfterMD toolkit.

This script demonstrates how to use the various modules of AfterMD
for MD trajectory analysis.
"""

import os
from pathlib import Path
from aftermd import (
    BatchProcessor, 
    PathManager, 
    PlotManager, 
    PBCProcessor, 
    TrajectoryAnalyzer, 
    StructureAnalyzer
)


def example_path_management():
    """Example of using PathManager."""
    print("=== Path Management Example ===")
    
    # Initialize path manager
    pm = PathManager("/path/to/project")
    
    # Create project structure
    pm.create_project_structure()
    
    # Get specific paths
    results_path = pm.get_path('results')
    plots_path = pm.get_path('plots')
    
    # Get timestamped output path
    timestamped_path = pm.get_timestamped_path('analysis', 'rmsd_analysis')
    
    print(f"Results path: {results_path}")
    print(f"Plots path: {plots_path}")
    print(f"Timestamped path: {timestamped_path}")
    
    # Print path structure
    print(pm)


def example_trajectory_analysis():
    """Example of trajectory analysis."""
    print("\n=== Trajectory Analysis Example ===")
    
    # Initialize analyzer (replace with actual file paths)
    topology = "protein.tpr"
    trajectory = "trajectory.xtc"
    
    if os.path.exists(topology) and os.path.exists(trajectory):
        analyzer = TrajectoryAnalyzer(topology, trajectory)
        
        # Calculate RMSD
        times, rmsd = analyzer.calculate_rmsd(
            selection="protein and name CA",
            output_file="rmsd_analysis.csv"
        )
        
        # Calculate radius of gyration
        times_rg, rg = analyzer.calculate_gyration_radius(
            selection="protein",
            output_file="radius_gyration.csv"
        )
        
        # Get trajectory info
        info = analyzer.get_trajectory_info()
        print(f"Trajectory info: {info}")
        
    else:
        print("Trajectory files not found - using placeholder")


def example_structure_analysis():
    """Example of structure analysis."""
    print("\n=== Structure Analysis Example ===")
    
    # Initialize analyzer (replace with actual PDB file)
    structure_file = "protein.pdb"
    
    if os.path.exists(structure_file):
        analyzer = StructureAnalyzer(structure_file)
        
        # Extract B-factors
        indices, bfactors = analyzer.extract_bfactors(
            selection="protein",
            output_file="bfactors.csv"
        )
        
        # Analyze B-factors by residue
        residue_stats = analyzer.analyze_bfactor_by_residue(
            output_file="bfactor_by_residue.csv"
        )
        
        # Calculate contact map
        contact_map = analyzer.calculate_contact_map(
            selection="name CA",
            cutoff=8.0,
            output_file="contact_map.txt"
        )
        
        # Get structure info
        info = analyzer.get_structure_info()
        print(f"Structure info: {info}")
        
    else:
        print("Structure file not found - using placeholder")


def example_plotting():
    """Example of plotting capabilities."""
    print("\n=== Plotting Example ===")
    
    # Initialize plot manager
    plotter = PlotManager()
    
    # Example XVG file plotting (replace with actual file)
    xvg_file = "rmsd.xvg"
    
    if os.path.exists(xvg_file):
        # Plot time series
        plotter.plot_time_series(
            xvg_file,
            title="RMSD Over Time",
            xlabel="Time (ps)",
            ylabel="RMSD (nm)",
            output_path="rmsd_plot.png"
        )
        
        # Plot RMSD specifically
        plotter.plot_rmsd(xvg_file, output_path="rmsd_analysis.png")
        
    else:
        print("XVG file not found - skipping plotting example")


def example_preprocessing():
    """Example of PBC processing workflow."""
    print("\n=== PBC Processing Example ===")
    
    # Initialize PBC processor
    pbc_processor = PBCProcessor()
    
    # Example files (replace with actual paths)
    trajectory = "raw_trajectory.xtc"
    topology = "system.tpr"
    
    if os.path.exists(trajectory) and os.path.exists(topology):
        # Remove PBC with trajectory downsampling
        processed_traj = pbc_processor.remove_pbc(
            trajectory=trajectory,
            topology=topology,
            output="processed_trajectory.xtc",
            dt=10.0  # Sample every 10 ps
        )
        
        # Comprehensive PBC processing with downsampling
        results = pbc_processor.comprehensive_pbc_process(
            trajectory=trajectory,
            topology=topology,
            output_dir="processed_data",
            dt=20.0  # Sample every 20 ps
        )
        
        print(f"PBC processing results: {results}")
        
    else:
        print("Input files not found - using placeholder")


def example_batch_processing():
    """Example of batch processing."""
    print("\n=== Batch Processing Example ===")
    
    # Initialize batch processor
    processor = BatchProcessor(max_workers=4)
    
    # Find trajectory files
    trajectory_files = processor.find_files(".", "*.xtc")
    print(f"Found {len(trajectory_files)} trajectory files")
    
    # Example processing function
    def analyze_trajectory(traj_file, output_dir="results"):
        """Example analysis function for batch processing."""
        try:
            # This would be replaced with actual analysis
            print(f"Processing: {traj_file}")
            return f"Processed {Path(traj_file).name}"
        except Exception as e:
            print(f"Error processing {traj_file}: {e}")
            return None
    
    if trajectory_files:
        # Batch analyze
        results = processor.batch_analyze(
            input_pattern="*.xtc",
            output_dir="batch_results",
            analysis_func=analyze_trajectory
        )
        
        print(f"Batch processing summary: {results}")
    else:
        print("No trajectory files found for batch processing")


def main():
    """Run all examples."""
    print("AfterMD Toolkit - Basic Usage Examples")
    print("=" * 50)
    
    example_path_management()
    example_trajectory_analysis()
    example_structure_analysis()
    example_plotting()
    example_preprocessing()
    example_batch_processing()
    
    print("\n" + "=" * 50)
    print("Examples completed!")
    print("Note: Replace placeholder file paths with actual data files")


if __name__ == "__main__":
    main()