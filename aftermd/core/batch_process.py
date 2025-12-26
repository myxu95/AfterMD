#!/usr/bin/env python3
"""
AfterMD Batch Processing - Simple Interface

A simple, direct interface for batch processing multiple MD simulation tasks.
Just provide a path containing task folders, and the system will automatically
discover and process all complete MD product files.

Usage:
    from aftermd.batch_process import process_md_tasks
    
    # Simple batch processing
    results = process_md_tasks("/path/to/simulations")
    
    # With options
    results = process_md_tasks(
        "/path/to/simulations",
        output_dir="/path/to/output",
        dt=10.0,
        max_workers=4
    )
"""

import logging
from pathlib import Path
from typing import Dict, Any, Optional
from ..utils.batch_processor import BatchProcessor

logger = logging.getLogger(__name__)


def process_md_tasks(simulations_path: str,
                    output_dir: Optional[str] = None,
                    dt: Optional[float] = None,
                    max_workers: Optional[int] = None,
                    gmx_executable: str = "gmx") -> Dict[str, Any]:
    """
    Process multiple MD simulation tasks with automatic file discovery.
    
    This function provides a simple interface for batch processing MD simulations.
    It automatically discovers all valid tasks in the provided directory and
    processes them using the AfterMD PBC processing pipeline with shortest
    chain detection.
    
    Args:
        simulations_path: Path containing multiple task folders
        output_dir: Output directory (default: simulations_path + "_processed")
        dt: Time interval for frame sampling in ps (optional)
        max_workers: Maximum number of parallel workers (optional)
        gmx_executable: GROMACS executable command (default: "gmx")
        
    Returns:
        Dictionary with processing results and summary
        
    Directory Structure Expected:
        simulations_path/
        ├── task1/
        │   ├── md.xtc
        │   └── md.tpr
        ├── task2/
        │   └── prod/
        │       ├── md.xtc
        │       └── md.tpr
        └── task3/
            ├── md.xtc
            └── md.tpr
            
    Example:
        >>> results = process_md_tasks("/data/my_simulations")
        >>> print(f"Processed {results['successful']}/{results['total_tasks']} tasks")
    """
    
    # Validate input path
    sims_path = Path(simulations_path)
    if not sims_path.exists():
        raise ValueError(f"Simulations path does not exist: {simulations_path}")
    
    if not sims_path.is_dir():
        raise ValueError(f"Simulations path is not a directory: {simulations_path}")
    
    # Set default output directory
    if output_dir is None:
        output_dir = str(sims_path.parent / f"{sims_path.name}_processed")
    
    logger.info(f"Starting batch MD processing")
    logger.info(f"Input directory: {simulations_path}")
    logger.info(f"Output directory: {output_dir}")
    if dt:
        logger.info(f"Frame sampling: {dt} ps")
    if max_workers:
        logger.info(f"Max workers: {max_workers}")
    
    # Initialize batch processor
    batch_processor = BatchProcessor(max_workers=max_workers)
    
    # Run batch PBC processing
    try:
        results = batch_processor.batch_pbc_processing(
            base_directory=simulations_path,
            output_base_dir=output_dir,
            dt=dt,
            gmx_executable=gmx_executable
        )
        
        # Enhanced logging
        total = results['total_tasks']
        successful = results['successful']
        failed = results['failed']
        
        logger.info(f"Batch processing completed!")
        logger.info(f"Summary: {successful}/{total} tasks processed successfully")
        
        if failed > 0:
            logger.warning(f"{failed} tasks failed - check individual task logs for details")
        
        # Log individual task results
        for result in results['results']:
            task_name = result.get('task_name', 'unknown')
            status = result.get('status', 'unknown')
            
            if status == 'success':
                output_file = result.get('processed', 'unknown')
                logger.info(f"✅ {task_name}: {output_file}")
            else:
                error = result.get('error', 'unknown error')
                logger.error(f"❌ {task_name}: {error}")
        
        return results
        
    except Exception as e:
        logger.error(f"Batch processing failed: {e}")
        raise


def discover_md_tasks(simulations_path: str) -> Dict[str, tuple]:
    """
    Discover all valid MD tasks in a directory without processing them.
    
    Args:
        simulations_path: Path containing multiple task folders
        
    Returns:
        Dictionary mapping task names to (trajectory, topology) file pairs
        
    Example:
        >>> tasks = discover_md_tasks("/data/my_simulations")
        >>> for task_name, (traj, topo) in tasks.items():
        ...     print(f"{task_name}: {traj}, {topo}")
    """
    
    sims_path = Path(simulations_path)
    if not sims_path.exists():
        raise ValueError(f"Simulations path does not exist: {simulations_path}")
    
    batch_processor = BatchProcessor()
    discovered_tasks = batch_processor.discover_batch_tasks(simulations_path)
    
    logger.info(f"Discovered {len(discovered_tasks)} valid MD tasks in {simulations_path}")
    
    return discovered_tasks


def check_task_status(simulations_path: str) -> Dict[str, Dict[str, Any]]:
    """
    Check the status of MD tasks in a directory (valid/invalid and reasons).
    
    Args:
        simulations_path: Path containing multiple task folders
        
    Returns:
        Dictionary with detailed status information for each task
        
    Example:
        >>> status = check_task_status("/data/my_simulations")
        >>> for task, info in status.items():
        ...     print(f"{task}: {info['status']} - {info['reason']}")
    """
    
    sims_path = Path(simulations_path)
    if not sims_path.exists():
        raise ValueError(f"Simulations path does not exist: {simulations_path}")
    
    batch_processor = BatchProcessor()
    task_status = {}
    
    # Check each subdirectory
    for item in sims_path.iterdir():
        if item.is_dir():
            task_name = item.name
            md_files = batch_processor.find_md_input_files(str(item))
            
            if md_files:
                trajectory, topology = md_files
                traj_path = Path(trajectory)
                
                location = "prod subfolder" if traj_path.parent.name == "prod" else "task directory"
                
                task_status[task_name] = {
                    'status': 'valid',
                    'reason': f'Complete MD files found in {location}',
                    'trajectory': trajectory,
                    'topology': topology,
                    'location': location
                }
            else:
                # Analyze what's missing
                direct_xtc = item / "md.xtc"
                direct_tpr = item / "md.tpr"
                prod_xtc = item / "prod" / "md.xtc"
                prod_tpr = item / "prod" / "md.tpr"
                
                missing_files = []
                if not direct_xtc.exists() and not prod_xtc.exists():
                    missing_files.append("md.xtc")
                if not direct_tpr.exists() and not prod_tpr.exists():
                    missing_files.append("md.tpr")
                
                if missing_files:
                    reason = f"Missing files: {', '.join(missing_files)}"
                else:
                    reason = "Files exist but not in expected locations"
                
                task_status[task_name] = {
                    'status': 'invalid',
                    'reason': reason,
                    'trajectory': None,
                    'topology': None,
                    'location': None
                }
    
    logger.info(f"Checked {len(task_status)} tasks in {simulations_path}")
    
    return task_status


# Convenience function for command line usage
def main():
    """Command line interface for batch processing."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="AfterMD Batch Processing - Process multiple MD simulation tasks",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python -m aftermd.batch_process /data/simulations
    python -m aftermd.batch_process /data/simulations --output /data/processed
    python -m aftermd.batch_process /data/simulations --dt 10.0 --workers 4
    python -m aftermd.batch_process /data/simulations --check-only
        """
    )
    
    parser.add_argument(
        "simulations_path",
        help="Path containing multiple task folders with MD files"
    )
    
    parser.add_argument(
        "--output", "-o",
        help="Output directory (default: simulations_path + '_processed')"
    )
    
    parser.add_argument(
        "--dt",
        type=float,
        help="Time interval for frame sampling in ps"
    )
    
    parser.add_argument(
        "--workers", "-w",
        type=int,
        help="Maximum number of parallel workers"
    )
    
    parser.add_argument(
        "--gmx",
        default="gmx",
        help="GROMACS executable command (default: gmx)"
    )
    
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Only check task status, don't process"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    try:
        if args.check_only:
            # Just check status
            print(f"Checking MD tasks in: {args.simulations_path}")
            print()
            
            status = check_task_status(args.simulations_path)
            
            valid_count = sum(1 for info in status.values() if info['status'] == 'valid')
            total_count = len(status)
            
            print(f"Task Status Summary: {valid_count}/{total_count} valid tasks")
            print("=" * 60)
            
            for task_name, info in status.items():
                status_icon = "✅" if info['status'] == 'valid' else "❌"
                print(f"{status_icon} {task_name}: {info['reason']}")
                
                if info['status'] == 'valid':
                    print(f"   📁 {info['trajectory']}")
                    print(f"   📁 {info['topology']}")
                print()
        
        else:
            # Run processing
            results = process_md_tasks(
                simulations_path=args.simulations_path,
                output_dir=args.output,
                dt=args.dt,
                max_workers=args.workers,
                gmx_executable=args.gmx
            )
            
            print()
            print("Batch Processing Results:")
            print("=" * 40)
            print(f"Total tasks: {results['total_tasks']}")
            print(f"Successful: {results['successful']}")
            print(f"Failed: {results['failed']}")
            print(f"Output directory: {results['output_directory']}")
            
    except Exception as e:
        logger.error(f"Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())