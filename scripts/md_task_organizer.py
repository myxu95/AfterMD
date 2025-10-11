#!/usr/bin/env python3
"""
MD Task Organizer Script

This script organizes MD simulation tasks by:
1. Finding incomplete tasks (missing md.xtc or md.tpr files)
2. Moving incomplete tasks to a separate folder
3. Comparing PDB IDs between incomplete and completed tasks
4. Removing duplicates (keeping completed versions)

Usage:
    python scripts/md_task_organizer.py /path/to/md/tasks --incomplete-dir /path/to/incomplete
    python scripts/md_task_organizer.py /path/to/md/tasks --dry-run
"""

import argparse
import shutil
import logging
from pathlib import Path
from typing import List, Dict, Set, Tuple
import re

logger = logging.getLogger(__name__)


def extract_pdb_id(task_name: str) -> str:
    """
    Extract PDB ID from task name.

    Args:
        task_name: Task directory name (e.g., "1abc_1", "2def_processed", "3ghi_2")

    Returns:
        PDB ID (first 4 characters in lowercase)
    """
    # Extract first 4 characters that could be a PDB ID
    match = re.match(r'^([a-zA-Z0-9]{4})', task_name)
    if match:
        return match.group(1).lower()
    else:
        # Fallback: use first 4 characters
        return task_name[:4].lower()


def check_md_files_exist(task_dir: Path) -> Tuple[bool, str]:
    """
    Check if MD files exist in task directory or prod subfolder.

    Args:
        task_dir: Task directory path

    Returns:
        Tuple of (files_exist, location)
    """
    # Check in task directory directly
    direct_xtc = task_dir / "md.xtc"
    direct_tpr = task_dir / "md.tpr"

    if direct_xtc.exists() and direct_tpr.exists():
        return True, "task_directory"

    # Check in prod subfolder
    prod_dir = task_dir / "prod"
    if prod_dir.exists():
        prod_xtc = prod_dir / "md.xtc"
        prod_tpr = prod_dir / "md.tpr"

        if prod_xtc.exists() and prod_tpr.exists():
            return True, "prod_subfolder"

    return False, "missing"


def find_incomplete_tasks(root_dir: Path) -> List[Dict[str, str]]:
    """
    Find all incomplete MD tasks in the root directory.

    Args:
        root_dir: Root directory containing MD tasks

    Returns:
        List of incomplete task information
    """
    incomplete_tasks = []

    for item in root_dir.iterdir():
        if not item.is_dir():
            continue

        task_name = item.name
        files_exist, location = check_md_files_exist(item)

        if not files_exist:
            pdb_id = extract_pdb_id(task_name)
            incomplete_tasks.append({
                "name": task_name,
                "path": str(item),
                "pdb_id": pdb_id
            })
            logger.debug(f"Found incomplete task: {task_name} (PDB: {pdb_id})")
        else:
            logger.debug(f"Task {task_name} is complete in {location}")

    return incomplete_tasks


def find_completed_tasks(root_dir: Path) -> Set[str]:
    """
    Find PDB IDs of all completed MD tasks.

    Args:
        root_dir: Root directory containing MD tasks

    Returns:
        Set of PDB IDs that have completed tasks
    """
    completed_pdb_ids = set()

    for item in root_dir.iterdir():
        if not item.is_dir():
            continue

        task_name = item.name
        files_exist, location = check_md_files_exist(item)

        if files_exist:
            pdb_id = extract_pdb_id(task_name)
            completed_pdb_ids.add(pdb_id)
            logger.debug(f"Completed task: {task_name} (PDB: {pdb_id}) in {location}")

    return completed_pdb_ids


def move_incomplete_tasks(incomplete_tasks: List[Dict[str, str]],
                         incomplete_dir: Path,
                         dry_run: bool = False) -> Dict[str, int]:
    """
    Move incomplete tasks to the incomplete directory.

    Args:
        incomplete_tasks: List of incomplete task information
        incomplete_dir: Directory to move incomplete tasks to
        dry_run: If True, only simulate the operation

    Returns:
        Dictionary with operation statistics
    """
    if not dry_run:
        incomplete_dir.mkdir(parents=True, exist_ok=True)

    stats = {"moved": 0, "failed": 0, "skipped": 0}

    for task in incomplete_tasks:
        source_path = Path(task["path"])
        dest_path = incomplete_dir / task["name"]

        try:
            if dest_path.exists():
                logger.warning(f"Destination already exists, skipping: {dest_path}")
                stats["skipped"] += 1
                continue

            if dry_run:
                logger.info(f"[DRY RUN] Would move: {source_path} -> {dest_path}")
            else:
                shutil.move(str(source_path), str(dest_path))
                logger.info(f"Moved: {source_path} -> {dest_path}")

            stats["moved"] += 1

        except Exception as e:
            logger.error(f"Failed to move {source_path}: {e}")
            stats["failed"] += 1

    return stats


def remove_duplicate_incomplete_tasks(incomplete_tasks: List[Dict[str, str]],
                                    completed_pdb_ids: Set[str],
                                    dry_run: bool = False) -> Dict[str, int]:
    """
    Remove incomplete tasks that have corresponding completed versions.

    Args:
        incomplete_tasks: List of incomplete task information
        completed_pdb_ids: Set of PDB IDs with completed tasks
        dry_run: If True, only simulate the operation

    Returns:
        Dictionary with operation statistics
    """
    stats = {"removed": 0, "failed": 0, "kept": 0}

    for task in incomplete_tasks:
        pdb_id = task["pdb_id"]
        task_path = Path(task["path"])

        if pdb_id in completed_pdb_ids:
            # This incomplete task has a completed version, remove it
            try:
                if dry_run:
                    logger.info(f"[DRY RUN] Would remove duplicate: {task_path} (PDB: {pdb_id})")
                else:
                    if task_path.exists():
                        shutil.rmtree(task_path)
                        logger.info(f"Removed duplicate: {task_path} (PDB: {pdb_id})")
                    else:
                        logger.warning(f"Path does not exist: {task_path}")

                stats["removed"] += 1

            except Exception as e:
                logger.error(f"Failed to remove {task_path}: {e}")
                stats["failed"] += 1
        else:
            # Keep this incomplete task (no completed version found)
            logger.info(f"Keeping unique incomplete task: {task['name']} (PDB: {pdb_id})")
            stats["kept"] += 1

    return stats


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="MD Task Organizer - Organize incomplete and duplicate MD tasks",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Organize tasks and move incomplete ones
    python scripts/md_task_organizer.py /data/md_tasks --incomplete-dir /data/incomplete_tasks

    # Remove duplicates only (keep completed versions)
    python scripts/md_task_organizer.py /data/md_tasks --remove-duplicates-only

    # Dry run to see what would happen
    python scripts/md_task_organizer.py /data/md_tasks --incomplete-dir /data/incomplete --dry-run

    # Verbose output
    python scripts/md_task_organizer.py /data/md_tasks --incomplete-dir /data/incomplete --verbose
        """
    )

    # Required arguments
    parser.add_argument(
        "task_directory",
        help="Directory containing MD task folders"
    )

    # Organization options
    parser.add_argument(
        "--incomplete-dir",
        help="Directory to move incomplete tasks to"
    )

    parser.add_argument(
        "--remove-duplicates-only",
        action="store_true",
        help="Only remove duplicates, don't move incomplete tasks"
    )

    # Operation modes
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without making changes"
    )

    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging"
    )

    args = parser.parse_args()

    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Validate arguments
    task_dir = Path(args.task_directory)
    if not task_dir.exists():
        logger.error(f"Task directory does not exist: {task_dir}")
        return 1

    if not task_dir.is_dir():
        logger.error(f"Task directory is not a directory: {task_dir}")
        return 1

    if not args.remove_duplicates_only and not args.incomplete_dir:
        logger.error("Either --incomplete-dir or --remove-duplicates-only must be specified")
        return 1

    print("🔍 MD Task Organizer")
    print("=" * 50)
    print(f"📂 Task directory: {task_dir}")

    if args.dry_run:
        print("🧪 DRY RUN MODE - No changes will be made")

    # Find incomplete and completed tasks
    logger.info("Scanning for incomplete tasks...")
    incomplete_tasks = find_incomplete_tasks(task_dir)

    logger.info("Scanning for completed tasks...")
    completed_pdb_ids = find_completed_tasks(task_dir)

    print(f"\n📊 Task Analysis:")
    print(f"❌ Incomplete tasks: {len(incomplete_tasks)}")
    print(f"✅ Completed PDB IDs: {len(completed_pdb_ids)}")

    # Find duplicates
    duplicate_count = 0
    unique_incomplete = 0

    for task in incomplete_tasks:
        if task["pdb_id"] in completed_pdb_ids:
            duplicate_count += 1
        else:
            unique_incomplete += 1

    print(f"🔄 Duplicate incomplete tasks: {duplicate_count}")
    print(f"🆕 Unique incomplete tasks: {unique_incomplete}")

    if not incomplete_tasks:
        print("\n✨ No incomplete tasks found. All tasks appear to be complete!")
        return 0

    # Move incomplete tasks (if requested)
    if args.incomplete_dir and not args.remove_duplicates_only:
        incomplete_dir = Path(args.incomplete_dir)
        print(f"\n📦 Moving incomplete tasks to: {incomplete_dir}")

        move_stats = move_incomplete_tasks(incomplete_tasks, incomplete_dir, args.dry_run)

        print(f"   ✅ Moved: {move_stats['moved']}")
        print(f"   ❌ Failed: {move_stats['failed']}")
        print(f"   ⏭️  Skipped: {move_stats['skipped']}")

    # Remove duplicates
    print(f"\n🗑️  Removing duplicate incomplete tasks...")

    duplicate_stats = remove_duplicate_incomplete_tasks(
        incomplete_tasks, completed_pdb_ids, args.dry_run
    )

    print(f"   🗑️  Removed: {duplicate_stats['removed']}")
    print(f"   ✅ Kept unique: {duplicate_stats['kept']}")
    print(f"   ❌ Failed: {duplicate_stats['failed']}")

    print(f"\n✨ Organization completed!")

    if args.dry_run:
        print("💡 Run without --dry-run to apply changes")

    return 0


if __name__ == "__main__":
    exit(main())