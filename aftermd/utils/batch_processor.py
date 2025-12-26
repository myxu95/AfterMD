import os
import glob
from pathlib import Path
from typing import List, Callable, Any, Dict, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import logging

logger = logging.getLogger(__name__)


# Global function for multiprocessing compatibility
def _process_md_task_global(task_info: Tuple[str, Tuple[str, str]], 
                           output_base_dir: str, dt: Optional[float], 
                           gmx_executable: str) -> Dict[str, Any]:
    """
    Global function to process a single MD task - required for pickle compatibility.
    """
    from ..preprocessing.pbc_processor import PBCProcessor
    
    task_name, (trajectory, topology) = task_info
    
    try:
        logger.info(f"Processing task: {task_name}")
        
        # Create task-specific output directory
        task_output_dir = Path(output_base_dir) / task_name
        task_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize PBC processor
        pbc_processor = PBCProcessor(gmx_executable)
        
        # Run comprehensive PBC processing
        result = pbc_processor.comprehensive_pbc_process(
            trajectory=trajectory,
            topology=topology,
            output_dir=str(task_output_dir),
            dt=dt
        )
        
        result["task_name"] = task_name
        result["status"] = "success"
        logger.info(f"Successfully processed task: {task_name}")
        return result
        
    except Exception as e:
        logger.error(f"Failed to process task {task_name}: {e}")
        return {
            "task_name": task_name,
            "status": "failed",
            "error": str(e)
        }


class BatchProcessor:
    """Batch processing utility for handling multiple files efficiently."""
    
    def __init__(self, max_workers: int = None):
        """
        Initialize BatchProcessor.
        
        Args:
            max_workers: Maximum number of worker processes/threads
        """
        self.max_workers = max_workers
        
    def find_files(self, directory: str, pattern: str) -> List[str]:
        """
        Find files matching pattern in directory.
        
        Args:
            directory: Directory path to search
            pattern: File pattern (e.g., "*.xtc", "*.pdb")
            
        Returns:
            List of matching file paths
        """
        search_path = os.path.join(directory, pattern)
        files = glob.glob(search_path, recursive=True)
        return sorted(files)
    
    def find_md_input_files(self, task_path: str) -> Optional[Tuple[str, str]]:
        """
        Find MD input files (md.xtc and md.tpr) using specific discovery rules.
        
        Discovery rules:
        1. Look directly in the task path
        2. Look in "prod" subfolder under task path
        
        Args:
            task_path: Path to the task directory
            
        Returns:
            Tuple of (trajectory_file, topology_file) if both found, None otherwise
        """
        task_dir = Path(task_path)
        
        if not task_dir.exists():
            logger.warning(f"Task directory does not exist: {task_path}")
            return None
        
        # Define search locations in priority order
        search_locations = [
            task_dir,              # Rule 1: directly in task path
            task_dir / "prod"      # Rule 2: in prod subfolder
        ]
        
        for location in search_locations:
            if not location.exists():
                continue
                
            xtc_file = location / "md.xtc"
            tpr_file = location / "md.tpr"
            
            if xtc_file.exists() and tpr_file.exists():
                logger.info(f"Found MD files in {location}: md.xtc, md.tpr")
                return str(xtc_file), str(tpr_file)
            else:
                missing = []
                if not xtc_file.exists():
                    missing.append("md.xtc")
                if not tpr_file.exists():
                    missing.append("md.tpr")
                logger.debug(f"Location {location}: missing {', '.join(missing)}")
        
        logger.warning(f"No complete MD input files (md.xtc + md.tpr) found in task: {task_path}")
        return None
    
    def discover_batch_tasks(self, base_directory: str) -> Dict[str, Tuple[str, str]]:
        """
        Discover all valid MD tasks in a base directory.
        
        Args:
            base_directory: Base directory containing multiple task folders
            
        Returns:
            Dictionary mapping task names to (trajectory, topology) file pairs
        """
        base_path = Path(base_directory)
        
        if not base_path.exists():
            logger.error(f"Base directory does not exist: {base_directory}")
            return {}
        
        discovered_tasks = {}
        
        # Look for subdirectories that might contain MD files
        for item in base_path.iterdir():
            if item.is_dir():
                task_name = item.name
                md_files = self.find_md_input_files(str(item))

                if md_files:
                    trajectory, topology = md_files
                    discovered_tasks[task_name] = (trajectory, topology)
                    logger.info(f"Task '{task_name}': found MD files")
                else:
                    logger.debug(f"Task '{task_name}': no MD files found")
        
        logger.info(f"Discovered {len(discovered_tasks)} valid MD tasks in {base_directory}")
        return discovered_tasks
    
    def process_files(self, 
                     files: List[str], 
                     process_func: Callable,
                     use_multiprocessing: bool = True,
                     **kwargs) -> List[Any]:
        """
        Process multiple files using parallel execution.
        
        Args:
            files: List of file paths to process
            process_func: Function to apply to each file
            use_multiprocessing: Use ProcessPoolExecutor if True, ThreadPoolExecutor if False
            **kwargs: Additional arguments to pass to process_func
            
        Returns:
            List of results from processing each file
        """
        executor_class = ProcessPoolExecutor if use_multiprocessing else ThreadPoolExecutor
        
        with executor_class(max_workers=self.max_workers) as executor:
            futures = [executor.submit(process_func, file_path, **kwargs) for file_path in files]
            results = []
            
            for future in futures:
                try:
                    result = future.result()
                    results.append(result)
                except Exception as e:
                    logger.error(f"Error processing file: {e}")
                    results.append(None)
                    
        return results
    
    def batch_analyze(self, 
                     input_pattern: str,
                     output_dir: str,
                     analysis_func: Callable,
                     **analysis_kwargs) -> Dict[str, Any]:
        """
        Perform batch analysis on files matching pattern.
        
        Args:
            input_pattern: Pattern to find input files (e.g., "/path/to/*.xtc")
            output_dir: Directory to save results
            analysis_func: Analysis function to apply
            **analysis_kwargs: Additional arguments for analysis function
            
        Returns:
            Dictionary with analysis results summary
        """
        files = glob.glob(input_pattern)
        if not files:
            raise ValueError(f"No files found matching pattern: {input_pattern}")
            
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        results = self.process_files(files, analysis_func, 
                                   output_dir=output_dir, **analysis_kwargs)
        
        successful = sum(1 for r in results if r is not None)
        total = len(files)
        
        summary = {
            "total_files": total,
            "successful": successful,
            "failed": total - successful,
            "results": results
        }
        
        logger.info(f"Batch analysis completed: {successful}/{total} files processed successfully")
        return summary
    
    def batch_pbc_processing(self, 
                           base_directory: str, 
                           output_base_dir: str,
                           dt: Optional[float] = None,
                           gmx_executable: str = "gmx") -> Dict[str, Any]:
        """
        Perform batch PBC processing with shortest chain detection on multiple MD tasks.
        
        Args:
            base_directory: Base directory containing task folders
            output_base_dir: Base directory for output files
            dt: Time interval for frame sampling in ps (optional)
            gmx_executable: GROMACS executable command
            
        Returns:
            Dictionary with batch processing results summary
        """
        from ..preprocessing.pbc_processor import PBCProcessor
        
        # Discover all valid MD tasks
        discovered_tasks = self.discover_batch_tasks(base_directory)
        
        if not discovered_tasks:
            logger.error("No valid MD tasks found for batch processing")
            return {"total_tasks": 0, "successful": 0, "failed": 0, "results": []}
        
        # Create output base directory
        Path(output_base_dir).mkdir(parents=True, exist_ok=True)
        
        # Prepare task list for parallel processing
        task_list = list(discovered_tasks.items())
        
        # Use custom parallel processing for tasks with global function
        executor_class = ProcessPoolExecutor
        with executor_class(max_workers=self.max_workers) as executor:
            futures = [
                executor.submit(_process_md_task_global, task_info, output_base_dir, dt, gmx_executable) 
                for task_info in task_list
            ]
            results = []
            
            for future in futures:
                try:
                    result = future.result()
                    results.append(result)
                except Exception as e:
                    logger.error(f"Error in batch processing: {e}")
                    results.append({"status": "failed", "error": str(e)})
        
        # Summarize results
        successful = sum(1 for r in results if r.get("status") == "success")
        failed = len(results) - successful
        
        summary = {
            "total_tasks": len(discovered_tasks),
            "successful": successful,
            "failed": failed,
            "results": results,
            "output_directory": output_base_dir
        }
        
        logger.info(f"Batch PBC processing completed: {successful}/{len(discovered_tasks)} tasks processed successfully")
        
        # Log task details
        for result in results:
            task_name = result.get("task_name", "unknown")
            status = result.get("status", "unknown")
            if status == "success":
                logger.info(f"✅ {task_name}: PBC processing completed")
            else:
                error = result.get("error", "unknown error")
                logger.error(f"❌ {task_name}: {error}")
        
        return summary