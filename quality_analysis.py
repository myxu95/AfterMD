#!/usr/bin/env python3
"""
MD Production Quality Analysis Launcher

This script provides a command-line interface for running comprehensive quality analysis
on MD simulation production data using the AfterMD quality analysis modules.

Usage:
    python quality_analysis.py --input /path/to/md/data --output /path/to/results
    python quality_analysis.py --batch /path/to/batch/dirs --output /path/to/results
    python quality_analysis.py --pdb-dir /path/to/pdb/files --output /path/to/results
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Optional
import json

from aftermd.analysis.quality import (
    MDCompletenessChecker,
    StructureValidator,
    BatchTracker,
    QualityReporter
)


def setup_logging(log_level: str = "INFO", log_file: Optional[str] = None):
    """Setup logging configuration."""
    level = getattr(logging, log_level.upper())
    
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )


def run_md_completeness_analysis(input_path: Path, 
                                min_traj_size: float,
                                min_sim_time: float) -> dict:
    """Run MD completeness analysis."""
    logger = logging.getLogger(__name__)
    logger.info(f"Starting MD completeness analysis on: {input_path}")
    
    checker = MDCompletenessChecker(
        min_trajectory_size_mb=min_traj_size,
        min_simulation_time_ps=min_sim_time
    )
    
    if input_path.is_dir():
        results = checker.batch_check(str(input_path))
    else:
        results = checker.check_single_md(str(input_path))
    
    logger.info("MD completeness analysis completed")
    return results


def run_structure_validation(pdb_files: List[Path], 
                           expected_chains: int,
                           pdb_id_only: bool = True) -> dict:
    """Run structure validation analysis."""
    logger = logging.getLogger(__name__)
    logger.info(f"Starting structure validation on {len(pdb_files)} PDB files")
    
    validator = StructureValidator(expected_chain_count=expected_chains)
    results = validator.batch_validate([str(pdb) for pdb in pdb_files], pdb_id_only=pdb_id_only)
    
    logger.info("Structure validation completed")
    return results


def run_batch_tracking(batch_path: Path) -> dict:
    """Run batch tracking analysis."""
    logger = logging.getLogger(__name__)
    logger.info(f"Starting batch tracking analysis on: {batch_path}")
    
    tracker = BatchTracker()
    results = tracker.analyze_batch(str(batch_path))
    
    logger.info("Batch tracking analysis completed")
    return results


def generate_quality_report(completeness_results: dict,
                          validation_results: dict,
                          batch_results: dict,
                          output_dir: Path) -> None:
    """Generate comprehensive quality report."""
    logger = logging.getLogger(__name__)
    logger.info("Generating comprehensive quality report")
    
    reporter = QualityReporter()
    report = reporter.generate_comprehensive_report(
        completeness_results,
        validation_results, 
        batch_results
    )
    
    # Save detailed report
    report_file = output_dir / "quality_analysis_report.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    # Generate summary report
    summary_file = output_dir / "quality_summary.txt"
    reporter.save_summary_report(report, str(summary_file))
    
    logger.info(f"Quality reports saved to: {output_dir}")


def is_pdb_id_filename(filename: str) -> bool:
    """
    Check if filename follows PDB ID naming convention.
    
    PDB ID format: 4-character alphanumeric code (e.g., 1ABC, 2XYZ, 3DEF)
    Valid patterns:
    - 1abc.pdb, 2XYZ.pdb (4 chars + .pdb)
    - 1abc_clean.pdb, 2XYZ_processed.pdb (4 chars + suffix + .pdb)
    """
    import re
    
    # Remove extension
    name_without_ext = filename.lower()
    if name_without_ext.endswith('.pdb'):
        name_without_ext = name_without_ext[:-4]
    elif name_without_ext.endswith('.ent'):
        name_without_ext = name_without_ext[:-4]
    else:
        return False
    
    # Check for PDB ID pattern: digit + 3 alphanumeric (standard PDB format)
    pdb_id_pattern = r'^[0-9][a-z0-9]{3}(_.*)?$'
    
    # Additional validation: ensure it's exactly 4 characters before any underscore
    base_part = name_without_ext.split('_')[0]
    if len(base_part) != 4:
        return False
        
    return bool(re.match(pdb_id_pattern, name_without_ext))


def find_pdb_files(directory: Path, pdb_id_only: bool = True, recursive: bool = True) -> List[Path]:
    """
    Find PDB files in directory.
    
    Args:
        directory: Directory to search
        pdb_id_only: If True, only return files with PDB ID naming convention
        recursive: If True, search recursively; if False, only search in root directory
    """
    pdb_files = []
    
    for pattern in ["*.pdb", "*.PDB"]:
        if recursive:
            pdb_files.extend(directory.rglob(pattern))
        else:
            pdb_files.extend(directory.glob(pattern))
    
    if pdb_id_only:
        # Filter to only include PDB ID named files
        filtered_files = []
        for pdb_file in pdb_files:
            if is_pdb_id_filename(pdb_file.name):
                filtered_files.append(pdb_file)
        return sorted(filtered_files)
    
    return sorted(pdb_files)


def main():
    parser = argparse.ArgumentParser(
        description="MD Production Quality Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Analyze single MD directory
    python quality_analysis.py --input /path/to/md/run --output ./results
    
    # Analyze multiple MD directories
    python quality_analysis.py --batch /path/to/batch/root --output ./results
    
    # Analyze PDB structures only (PDB ID named files by default)
    python quality_analysis.py --pdb-dir /path/to/pdbs --output ./results
    
    # Include all PDB files (not just PDB ID named)
    python quality_analysis.py --pdb-dir /path/to/pdbs --include-all-pdb --output ./results
    
    # Custom parameters
    python quality_analysis.py --input /path/to/md \\
        --min-traj-size 5.0 --min-sim-time 10000 \\
        --expected-chains 3 --output ./results
        """
    )
    
    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--input", "-i", type=Path,
                           help="Path to single MD simulation directory")
    input_group.add_argument("--batch", "-b", type=Path,
                           help="Path to directory containing multiple MD runs")
    input_group.add_argument("--pdb-dir", "-p", type=Path,
                           help="Path to directory containing PDB files")
    
    # Output
    parser.add_argument("--output", "-o", type=Path, required=True,
                       help="Output directory for analysis results")
    
    # MD completeness parameters
    parser.add_argument("--min-traj-size", type=float, default=1.0,
                       help="Minimum trajectory file size in MB (default: 1.0)")
    parser.add_argument("--min-sim-time", type=float, default=5000.0,
                       help="Minimum simulation time in ps (default: 5000.0)")
    
    # Structure validation parameters
    parser.add_argument("--expected-chains", type=int, default=5,
                       help="Expected number of protein chains (default: 5)")
    
    # Logging options
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                       default="INFO", help="Logging level (default: INFO)")
    parser.add_argument("--log-file", type=Path,
                       help="Log file path (default: output to console)")
    
    # Analysis options
    parser.add_argument("--skip-completeness", action="store_true",
                       help="Skip MD completeness analysis")
    parser.add_argument("--skip-structure", action="store_true",
                       help="Skip structure validation")
    parser.add_argument("--skip-batch", action="store_true",
                       help="Skip batch tracking")
    
    # PDB file filtering options
    parser.add_argument("--include-all-pdb", action="store_true",
                       help="Include all PDB files, not just PDB ID named files")
    
    args = parser.parse_args()
    
    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)
    
    # Setup logging
    log_file = str(args.log_file) if args.log_file else None
    setup_logging(args.log_level, log_file)
    
    logger = logging.getLogger(__name__)
    logger.info("Starting MD Production Quality Analysis")
    logger.info(f"Output directory: {args.output}")
    
    # Initialize results
    completeness_results = {}
    validation_results = {}
    batch_results = {}
    
    try:
        # Determine analysis targets
        if args.input:
            analysis_path = args.input
            logger.info(f"Single MD analysis mode: {analysis_path}")
        elif args.batch:
            analysis_path = args.batch
            logger.info(f"Batch analysis mode: {analysis_path}")
        else:  # args.pdb_dir
            analysis_path = args.pdb_dir
            logger.info(f"PDB-only analysis mode: {analysis_path}")
        
        # Run MD completeness analysis
        if not args.skip_completeness and args.pdb_dir is None:
            completeness_results = run_md_completeness_analysis(
                analysis_path, args.min_traj_size, args.min_sim_time
            )
        
        # Run structure validation
        if not args.skip_structure:
            pdb_id_only = not args.include_all_pdb
            
            if args.pdb_dir:
                pdb_files = find_pdb_files(args.pdb_dir, pdb_id_only=pdb_id_only)
            else:
                # Find PDB files in MD Product directories (recursive search for MD Products)
                pdb_files = find_pdb_files(analysis_path, pdb_id_only=pdb_id_only, recursive=True)
            
            if pdb_files:
                if pdb_id_only:
                    logger.info(f"Found {len(pdb_files)} PDB files with PDB ID naming convention")
                else:
                    logger.info(f"Found {len(pdb_files)} PDB files (including all formats)")
                
                validation_results = run_structure_validation(
                    pdb_files, args.expected_chains, pdb_id_only=pdb_id_only
                )
            else:
                if pdb_id_only:
                    logger.warning("No PDB files with PDB ID naming convention found for structure validation")
                    logger.info("Use --include-all-pdb to include all PDB files")
                else:
                    logger.warning("No PDB files found for structure validation")
        
        # Run batch tracking
        if not args.skip_batch and args.pdb_dir is None:
            batch_results = run_batch_tracking(analysis_path)
        
        # Generate comprehensive report
        generate_quality_report(
            completeness_results,
            validation_results,
            batch_results,
            args.output
        )
        
        logger.info("Quality analysis completed successfully")
        
        # Print summary
        print("\n" + "="*60)
        print("MD PRODUCTION QUALITY ANALYSIS SUMMARY")
        print("="*60)
        
        if completeness_results:
            total_md = len(completeness_results.get('results', []))
            completed = sum(1 for r in completeness_results.get('results', []) 
                          if r.get('status') == 'complete')
            print(f"MD Completeness: {completed}/{total_md} simulations complete")
        
        if validation_results:
            total_pdb = len(validation_results.get('results', []))
            valid = sum(1 for r in validation_results.get('results', [])
                       if r.get('status') == 'valid')
            print(f"Structure Validation: {valid}/{total_pdb} structures valid")
        
        if batch_results:
            unique_pdbs = batch_results.get('summary', {}).get('unique_pdbs', 0)
            print(f"Batch Tracking: {unique_pdbs} unique PDB entries processed")
        
        print(f"\nDetailed results saved to: {args.output}")
        print("="*60)
        
    except Exception as e:
        logger.error(f"Quality analysis failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()