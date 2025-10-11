#!/usr/bin/env python3
"""
PBC Processing Startup Script

This script provides a command-line interface for PBC (Periodic Boundary Conditions)
processing using the AfterMD toolkit. It supports both single trajectory processing
and batch processing of multiple trajectory files.

Usage:
    python scripts/pbc_process.py -f trajectory.xtc -s topology.tpr -o output/
    python scripts/pbc_process.py -d input_dir/ -o output_dir/ --batch
"""

import argparse
import sys
import logging
from pathlib import Path
from typing import Optional, List

try:
    from aftermd.preprocessing import PBCProcessor
    from aftermd.utils import BatchProcessor, PathManager
except ImportError as e:
    print(f"Error importing AfterMD modules: {e}")
    print("Please ensure AfterMD is properly installed: pip install -e .")
    sys.exit(1)


def setup_logging(verbose: bool = False):
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def validate_input_files(trajectory: str, topology: str) -> bool:
    """Validate input files exist and are readable."""
    traj_path = Path(trajectory)
    topo_path = Path(topology)

    if not traj_path.exists():
        print(f"Error: Trajectory file not found: {trajectory}")
        return False

    if not topo_path.exists():
        print(f"Error: Topology file not found: {topology}")
        return False

    if not traj_path.is_file():
        print(f"Error: Trajectory path is not a file: {trajectory}")
        return False

    if not topo_path.is_file():
        print(f"Error: Topology path is not a file: {topology}")
        return False

    # Check file extensions
    valid_traj_exts = ['.xtc', '.trr']
    valid_topo_exts = ['.tpr', '.gro']

    if traj_path.suffix not in valid_traj_exts:
        print(f"Warning: Trajectory file extension {traj_path.suffix} not in {valid_traj_exts}")

    if topo_path.suffix not in valid_topo_exts:
        print(f"Warning: Topology file extension {topo_path.suffix} not in {valid_topo_exts}")

    return True


def process_single_trajectory(args) -> bool:
    """Process a single trajectory file."""
    print(f"Processing single trajectory: {args.trajectory}")
    print(f"Topology: {args.topology}")
    print(f"Output directory: {args.output}")

    # Validate input files
    if not validate_input_files(args.trajectory, args.topology):
        return False

    try:
        # Initialize PBC processor
        pbc_processor = PBCProcessor(
            gmx_executable=args.gmx,
            keep_temp_files=args.keep_temp
        )

        # Process trajectory
        results = pbc_processor.comprehensive_pbc_process(
            trajectory=args.trajectory,
            topology=args.topology,
            output_dir=args.output,
            center_group=args.center_group,
            fit_group=args.fit_group,
            dt=args.dt
        )

        print("PBC processing completed successfully!")
        print(f"Results: {results}")
        return True

    except Exception as e:
        print(f"Error during PBC processing: {e}")
        logging.exception("Detailed error information:")
        return False


def process_batch_trajectories(args) -> bool:
    """Process multiple trajectory files in batch mode."""
    print(f"Processing batch trajectories from: {args.input_dir}")
    print(f"Output directory: {args.output}")

    input_path = Path(args.input_dir)
    if not input_path.exists() or not input_path.is_dir():
        print(f"Error: Input directory not found: {args.input_dir}")
        return False

    try:
        # Initialize processors
        batch_processor = BatchProcessor(max_workers=args.workers)
        pbc_processor = PBCProcessor(
            gmx_executable=args.gmx,
            keep_temp_files=args.keep_temp
        )

        # Find trajectory files
        trajectory_files = batch_processor.find_files(
            str(input_path),
            args.pattern
        )

        if not trajectory_files:
            print(f"No trajectory files found matching pattern: {args.pattern}")
            return False

        print(f"Found {len(trajectory_files)} trajectory files")

        # Define processing function
        def process_trajectory(traj_file: str) -> Optional[dict]:
            """Process a single trajectory in batch mode."""
            traj_path = Path(traj_file)

            # Look for topology file in the same directory
            topo_candidates = [
                traj_path.parent / "md.tpr",
                traj_path.parent / "npt.tpr",
                traj_path.parent / "system.tpr",
                traj_path.parent / f"{traj_path.stem}.tpr"
            ]

            topology = None
            for topo_file in topo_candidates:
                if topo_file.exists():
                    topology = str(topo_file)
                    break

            if not topology:
                print(f"Warning: No topology file found for {traj_file}")
                return None

            # Create output directory for this trajectory
            output_subdir = Path(args.output) / f"{traj_path.stem}_processed"

            try:
                results = pbc_processor.comprehensive_pbc_process(
                    trajectory=traj_file,
                    topology=topology,
                    output_dir=str(output_subdir),
                    center_group=args.center_group,
                    fit_group=args.fit_group,
                    dt=args.dt
                )

                return {
                    'input': traj_file,
                    'output': results['processed'],
                    'status': 'success'
                }

            except Exception as e:
                print(f"Error processing {traj_file}: {e}")
                return {
                    'input': traj_file,
                    'error': str(e),
                    'status': 'failed'
                }

        # Execute batch processing
        print("Starting batch processing...")
        results = []

        for traj_file in trajectory_files:
            result = process_trajectory(traj_file)
            if result:
                results.append(result)

        # Summary
        successful = [r for r in results if r['status'] == 'success']
        failed = [r for r in results if r['status'] == 'failed']

        print(f"\nBatch processing completed!")
        print(f"Successful: {len(successful)}")
        print(f"Failed: {len(failed)}")

        if failed:
            print("\nFailed files:")
            for fail in failed:
                print(f"  {fail['input']}: {fail['error']}")

        return len(failed) == 0

    except Exception as e:
        print(f"Error during batch processing: {e}")
        logging.exception("Detailed error information:")
        return False


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="PBC Processing Script for MD Trajectories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process single trajectory
  python scripts/pbc_process.py -f md.xtc -s md.tpr -o processed/

  # Process with custom groups and downsampling
  python scripts/pbc_process.py -f md.xtc -s md.tpr -o processed/ --dt 10.0 --center-group "Protein" --fit-group "Backbone"

  # Batch process multiple trajectories
  python scripts/pbc_process.py -d input_dir/ -o output_dir/ --batch --pattern "*.xtc"

  # Keep temporary files for debugging
  python scripts/pbc_process.py -f md.xtc -s md.tpr -o processed/ --keep-temp --verbose
        """
    )

    # Input/output arguments
    parser.add_argument('-f', '--trajectory',
                        help='Input trajectory file (.xtc or .trr)')
    parser.add_argument('-s', '--topology',
                        help='Topology file (.tpr or .gro)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory')

    # Batch processing arguments
    parser.add_argument('-d', '--input-dir',
                        help='Input directory for batch processing')
    parser.add_argument('--batch', action='store_true',
                        help='Enable batch processing mode')
    parser.add_argument('--pattern', default='*.xtc',
                        help='File pattern for batch processing (default: *.xtc)')
    parser.add_argument('--workers', type=int, default=4,
                        help='Number of parallel workers for batch processing (default: 4)')

    # Processing options
    parser.add_argument('--center-group',
                        help='Group for centering (default: auto-select)')
    parser.add_argument('--fit-group',
                        help='Group for fitting (default: auto-select)')
    parser.add_argument('--dt', type=float,
                        help='Time interval for frame sampling in ps (default: no downsampling)')

    # GROMACS options
    parser.add_argument('--gmx', default='gmx',
                        help='GROMACS executable command (default: gmx)')

    # Debug options
    parser.add_argument('--keep-temp', action='store_true',
                        help='Keep temporary files for debugging')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose logging')

    args = parser.parse_args()

    # Setup logging
    setup_logging(args.verbose)

    # Validate arguments
    if args.batch:
        if not args.input_dir:
            print("Error: --input-dir required for batch processing")
            sys.exit(1)
    else:
        if not args.trajectory or not args.topology:
            print("Error: -f/--trajectory and -s/--topology required for single file processing")
            sys.exit(1)

    # Create output directory
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)

    print("AfterMD PBC Processing Script")
    print("=" * 40)

    # Execute processing
    if args.batch:
        success = process_batch_trajectories(args)
    else:
        success = process_single_trajectory(args)

    if success:
        print("\nProcessing completed successfully!")
        sys.exit(0)
    else:
        print("\nProcessing failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()