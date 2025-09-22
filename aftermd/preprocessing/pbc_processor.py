import subprocess
import os
from pathlib import Path
from typing import List, Optional, Dict, Any
import logging
from ..utils.group_selector import GroupSelector

logger = logging.getLogger(__name__)


class PBCProcessor:
    """PBC (Periodic Boundary Conditions) processing for MD trajectories."""
    
    def __init__(self, gmx_executable: str = "gmx"):
        """
        Initialize PBCProcessor.
        
        Args:
            gmx_executable: GROMACS executable command
        """
        self.gmx = gmx_executable
        self.group_selector = None  # Will be initialized when topology is provided
        self._check_gromacs()
    
    def _check_gromacs(self):
        """Check if GROMACS is available."""
        try:
            result = subprocess.run([self.gmx, "--version"], 
                                  capture_output=True, text=True, check=True)
            logger.info("GROMACS found and accessible")
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("GROMACS not found or not accessible")
    
    def _run_gmx_command(self, 
                        command: List[str], 
                        stdin_input: Optional[str] = None,
                        cwd: Optional[str] = None) -> subprocess.CompletedProcess:
        """
        Run GROMACS command with error handling.
        
        Args:
            command: GROMACS command list
            stdin_input: Input to pass to stdin
            cwd: Working directory
            
        Returns:
            CompletedProcess result
        """
        full_command = [self.gmx] + command
        logger.info(f"Running: {' '.join(full_command)}")
        
        try:
            result = subprocess.run(
                full_command,
                input=stdin_input,
                text=True,
                capture_output=True,
                check=True,
                cwd=cwd
            )
            return result
        except subprocess.CalledProcessError as e:
            logger.error(f"GROMACS command failed: {e}")
            logger.error(f"stdout: {e.stdout}")
            logger.error(f"stderr: {e.stderr}")
            raise
    
    def remove_pbc(self, 
                   trajectory: str,
                   topology: str,
                   output: str,
                   center_group: Optional[str] = None,
                   output_group: Optional[str] = None,
                   fit_group: Optional[str] = None,
                   dt: Optional[float] = None) -> str:
        """
        Remove PBC using three-step process: center -> pbc whole -> fit rot+trans.
        
        Args:
            trajectory: Input trajectory file (.xtc or .trr)
            topology: Topology file (.tpr or .gro)
            output: Output trajectory file
            center_group: Group to center (default: Protein)
            output_group: Group to output (default: System)
            fit_group: Group for fitting rot+trans (default: Protein)
            dt: Time interval for frame sampling in ps (default: no downsampling)
            
        Returns:
            Path to output file
        """
        # Initialize group selector for this topology
        if self.group_selector is None:
            self.group_selector = GroupSelector(topology, self.gmx)
        
        # Auto-select groups if not provided
        if center_group is None:
            center_group = self.group_selector.select_group("center_protein")
        if output_group is None:
            output_group = self.group_selector.select_group("output_system")
        if fit_group is None:
            fit_group = self.group_selector.select_group("fit_backbone")
        
        # Temporary files
        temp_centered = str(Path(output).with_suffix('.temp_centered.xtc'))
        temp_whole = str(Path(output).with_suffix('.temp_whole.xtc'))
        
        try:
            # Step 1: Center trajectory with PBC atom
            if dt is not None:
                logger.info(f"Step 1: Centering trajectory (center pbc atom) with dt={dt} ps")
            else:
                logger.info("Step 1: Centering trajectory (center pbc atom)")
            
            center_cmd = [
                "trjconv",
                "-f", trajectory,
                "-s", topology,
                "-o", temp_centered,
                "-center",
                "-pbc", "atom"
            ]
            
            # Add index file if shortest chain is available
            if self.group_selector and self.group_selector.has_shortest_chain_group():
                index_file = self.group_selector.get_shortest_chain_index_file()
                if index_file:
                    center_cmd.extend(["-n", index_file])
                    logger.info(f"Using custom index file for centering: {index_file}")
            
            # Add dt option if specified
            if dt is not None:
                center_cmd.extend(["-dt", str(dt)])
            
            center_input = f"{center_group}\n{output_group}\n"
            self._run_gmx_command(center_cmd, stdin_input=center_input)
            logger.info(f"Centered trajectory using group: {center_group}")
            
            # Step 2: Apply PBC whole
            logger.info("Step 2: Applying PBC whole")
            whole_cmd = [
                "trjconv",
                "-f", temp_centered,
                "-s", topology,
                "-o", temp_whole,
                "-pbc", "whole"
            ]
            
            # Note: dt not needed for step 2 since we're processing the already sampled trajectory
            
            whole_input = f"{output_group}\n"
            self._run_gmx_command(whole_cmd, stdin_input=whole_input)
            logger.info("Applied PBC whole")
            
            # Step 3: Fit rot+trans
            logger.info("Step 3: Fitting rotational and translational motion")
            fit_cmd = [
                "trjconv",
                "-f", temp_whole,
                "-s", topology,
                "-o", output,
                "-fit", "rot+trans"
            ]
            
            # Note: dt not needed for step 3 since we're processing the already sampled trajectory
            
            fit_input = f"{fit_group}\n{output_group}\n"
            self._run_gmx_command(fit_cmd, stdin_input=fit_input)
            logger.info(f"Applied rot+trans fitting using group: {fit_group}")
            
            logger.info(f"PBC processing completed: {output}")
            logger.info(f"Used groups - Center: {center_group}, Output: {output_group}, Fit: {fit_group}")
            return output
            
        finally:
            # Clean up temporary files
            temp_files = [temp_centered, temp_whole]
            
            # Add shortest chain index file to cleanup list if it exists
            if self.group_selector and self.group_selector.has_shortest_chain_group():
                index_file = self.group_selector.get_shortest_chain_index_file()
                if index_file:
                    temp_files.append(index_file)
            
            for temp_file in temp_files:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
                    logger.debug(f"Cleaned up temporary file: {temp_file}")
    
    def extract_frames(self,
                      trajectory: str,
                      topology: str,
                      output: str,
                      start_time: Optional[float] = None,
                      end_time: Optional[float] = None,
                      step: Optional[int] = None) -> str:
        """
        Extract specific frames from trajectory.
        
        Args:
            trajectory: Input trajectory file
            topology: Topology file
            output: Output trajectory file
            start_time: Start time in ps
            end_time: End time in ps
            step: Frame step interval
            
        Returns:
            Path to output file
        """
        cmd = [
            "trjconv",
            "-f", trajectory,
            "-s", topology,
            "-o", output
        ]
        
        if start_time is not None:
            cmd.extend(["-b", str(start_time)])
        if end_time is not None:
            cmd.extend(["-e", str(end_time)])
        if step is not None:
            cmd.extend(["-skip", str(step)])
        
        stdin_input = "0\n"
        self._run_gmx_command(cmd, stdin_input=stdin_input)
        
        logger.info(f"Extracted frames saved to: {output}")
        return output
    
    def comprehensive_pbc_process(self,
                                trajectory: str,
                                topology: str,
                                output_dir: str,
                                center_group: Optional[str] = None,
                                fit_group: Optional[str] = None,
                                dt: Optional[float] = None) -> Dict[str, str]:
        """
        Comprehensive PBC processing workflow using the new three-step process.
        
        Args:
            trajectory: Input trajectory file
            topology: Topology file
            output_dir: Output directory
            center_group: Group for centering (default: auto-select)
            fit_group: Group for fitting (default: auto-select)
            dt: Time interval for frame sampling in ps (default: no downsampling)
            
        Returns:
            Dictionary with paths to output files
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        base_name = Path(trajectory).stem
        
        # Initialize group selector
        if self.group_selector is None:
            self.group_selector = GroupSelector(topology, self.gmx)
        
        # Auto-select groups if not provided
        if center_group is None:
            center_group = self.group_selector.select_group("center_protein")
        if fit_group is None:
            fit_group = self.group_selector.select_group("fit_backbone")
        
        # Use the new three-step PBC removal process
        processed_file = output_path / f"{base_name}_processed.xtc"
        
        logger.info(f"Starting comprehensive PBC processing for {trajectory}")
        if dt is not None:
            logger.info(f"Using dt={dt} ps for trajectory downsampling")
        logger.info(f"Using groups - Center: {center_group}, Fit: {fit_group}")
        
        final_trajectory = self.remove_pbc(
            trajectory=trajectory,
            topology=topology,
            output=str(processed_file),
            center_group=center_group,
            fit_group=fit_group,
            dt=dt
        )
        
        results = {
            "processed": final_trajectory,
            "center_group": center_group,
            "fit_group": fit_group,
            "output_directory": str(output_path)
        }
        
        logger.info(f"Comprehensive PBC processing completed in: {output_dir}")
        return results