#!/bin/bash
#
# AfterMD Project Cleanup Script
# Removes redundant files and directories to keep the project clean
#

echo "=== AfterMD Project Cleanup ==="
echo "This script will remove redundant files and directories."
echo

# Function to safely remove files/directories
safe_remove() {
    local target="$1"
    local description="$2"
    
    if [ -e "$target" ]; then
        echo "Removing: $target ($description)"
        rm -rf "$target"
    else
        echo "Already removed: $target"
    fi
}

# Confirm before proceeding
read -p "Continue with cleanup? (y/N): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cleanup cancelled."
    exit 0
fi

echo
echo "Starting cleanup..."

# 1. Remove redundant documentation files
echo "--- Cleaning redundant documentation ---"
safe_remove "CLEAN_PROJECT_STRUCTURE.md" "redundant project structure doc"
safe_remove "FEATURES_SUMMARY.md" "redundant features summary"
safe_remove "INSTALLATION.md" "redundant installation doc"
safe_remove "PROJECT_OVERVIEW.md" "redundant project overview"
safe_remove "PROJECT_STATISTICS.md" "redundant statistics doc"
safe_remove "optimized_shortest_chain_guide.md" "outdated guide"
safe_remove "shortest_chain_detection_guide.md" "outdated guide"

# 2. Remove demo and test scripts
echo "--- Cleaning demo and test scripts ---"
safe_remove "demo_export_failed.py" "demo script"
safe_remove "demo_pdb_filtering.py" "demo script"
safe_remove "test_pdb_filtering.py" "test script"
safe_remove "test_structure_validator.py" "test script"

# 3. Remove outdated scripts
echo "--- Cleaning outdated scripts ---"
safe_remove "quick_quality_check.py" "replaced by md_quality_check.py"
safe_remove "pbc_process.sh" "standalone PBC script"
safe_remove "quick_slurm.sh" "simple SLURM script"
safe_remove "update2_28.sh" "old update script"

# 4. Remove empty/redundant directories
echo "--- Cleaning empty directories ---"
safe_remove "quality_reports" "empty reports directory"
safe_remove "docs" "redundant docs directory"

# 5. Clean Python cache files
echo "--- Cleaning Python cache ---"
find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
find . -name "*.pyc" -delete 2>/dev/null || true
find . -name "*.pyo" -delete 2>/dev/null || true

# 6. Clean temporary files
echo "--- Cleaning temporary files ---"
safe_remove "*.tmp" "temporary files"
safe_remove "*.log" "log files"
safe_remove ".DS_Store" "macOS system files"

echo
echo "=== Cleanup Summary ==="
echo "The following core files are preserved:"
echo "  • CLAUDE.md - Project configuration"
echo "  • README.md - Main project documentation"
echo "  • aftermd/ - Core Python package"
echo "  • examples/ - Usage examples"
echo "  • environment.yml - Conda environment"
echo "  • requirements.txt - Python dependencies"
echo "  • setup.py - Package setup"
echo
echo "  • md_quality_check.py - Quality checking tool"
echo "  • md_workflow.py - Complete workflow script"
echo "  • generate_slurm.py - SLURM script generator"
echo "  • quality_analysis.py - Quality analysis tool"
echo "  • export_failed_mds.py - Failed MD export tool"
echo "  • quality_aware_slurm_guide.md - Usage guide"
echo "  • slurm_template.sh - SLURM template"
echo
echo "✅ Project cleanup completed!"
echo "The project structure is now clean and focused."