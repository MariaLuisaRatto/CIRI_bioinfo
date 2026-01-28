#!/usr/bin/env python3

import os
import sys
import re
import subprocess
import pathlib
import logging
from typing import Dict, List, Any, Optional, Union
from dataclasses import dataclass

# Configure logging
logger = logging.getLogger(__name__)

@dataclass
class Result:
  status: str
  output_dir: str
  message: str = ""

def validate_path(path: str) -> str:
  """Validate and normalize a file path."""
  if not path:
    raise ValueError("Path cannot be empty")
  return os.path.abspath(os.path.expanduser(path))

def is_running_in_docker() -> bool:
  """Check if we're running inside a Docker container."""
  return os.path.exists('/.dockerenv')

def run_docker(image: str, volumes: Dict[str, str], env: Dict[str, str], args: List[str]) -> str:
  """Run a Docker container with specified parameters."""
  cmd = ['docker', 'run', '--rm']
  
  for src, dst in volumes.items():
    cmd.extend(['-v', f"{src}:{dst}"])
  
  for key, val in env.items():
    cmd.extend(['-e', f"{key}={val}"])
  
  cmd.append(image)
  cmd.extend(args)
  
  logger.info(f"Running Docker command: {' '.join(cmd)}")
  result = subprocess.run(cmd, capture_output=True, text=True, check=False)
  
  if result.returncode != 0:
    logger.error(f"Docker execution failed: {result.stderr}")
    raise RuntimeError(f"Docker execution failed: {result.stderr}")
  
  return result.stdout

def CIRI_assign(target_directory: str, script_directory: str, guide_strategy: str, matrix_filename: str, threshold_a: str, threshold_i: str) -> Result:
  """
  CRISPR QC and assignment analysis.
  
  Parameters:
      target_directory: Working directory containing the h5 matrix and guides.csv
      script_directory: Directory containing the CIRI_unified.R script
      guide_strategy: Variable guides logic: 1=Single, 2=Dual (allowed values: 1, 2)
      matrix_filename: Filename of the HDF5 feature matrix
      threshold_a: Manual threshold for CRISPRa
      threshold_i: Manual threshold for CRISPRi
  
  Returns:
      Result: Results of the operation
  """
  # Parameter validation
  if not isinstance(target_directory, str):
    raise TypeError(f"target_directory must be a string, got {type(target_directory).__name__}")
  target_directory_path = validate_path(target_directory)
  if not isinstance(script_directory, str):
    raise TypeError(f"script_directory must be a string, got {type(script_directory).__name__}")
  script_directory_path = validate_path(script_directory)
  guide_strategy_valid_values = ["1", "2"]
  if not isinstance(guide_strategy, str):
    raise TypeError(f"guide_strategy must be a string, got {type(guide_strategy).__name__}")
  if guide_strategy not in guide_strategy_valid_values:
    raise ValueError(f"guide_strategy must be one of {guide_strategy_valid_values}")
  if not isinstance(matrix_filename, str):
    raise TypeError(f"matrix_filename must be a string, got {type(matrix_filename).__name__}")
  if not isinstance(threshold_a, str):
    raise TypeError(f"threshold_a must be a string, got {type(threshold_a).__name__}")
  if not isinstance(threshold_i, str):
    raise TypeError(f"threshold_i must be a string, got {type(threshold_i).__name__}")
  
  # Directory existence checks
  if not is_running_in_docker():
    if not os.path.isdir(target_directory_path):
      raise NotADirectoryError(f"Directory {target_directory_path} does not exist")
  if not is_running_in_docker():
    if not os.path.isdir(script_directory_path):
      raise NotADirectoryError(f"Directory {script_directory_path} does not exist")
  
  # Process file paths for Docker volume mounting
  # Process target_directory for Docker
  target_directory_abspath = os.path.abspath(target_directory_path if 'target_directory_path' in locals() else target_directory)
  target_directory_dir = os.path.dirname(target_directory_abspath)
  target_directory_filename = os.path.basename(target_directory)
  # Process script_directory for Docker
  script_directory_abspath = os.path.abspath(script_directory_path if 'script_directory_path' in locals() else script_directory)
  script_directory_dir = os.path.dirname(script_directory_abspath)
  script_directory_filename = os.path.basename(script_directory)
  
  # Main volume mount point
  main_mount_dir = target_directory_dir
  
  # Execute Docker container with error handling
  try:
    # Prepare Docker volumes
    volumes = {}
    volumes[target_directory_dir] = "/mnt/analysis"
    volumes[script_directory_dir] = "/opt/scripts"
    
    # Prepare environment variables
    env_vars = {}
    
    # Prepare Docker arguments
    docker_args = []
    docker_args.append("/opt/scripts/CIRI_unified.R")
    docker_args.append("/mnt/analysis")
    docker_args.append(str(guide_strategy))
    docker_args.append(str(matrix_filename))
    docker_args.append(str(threshold_a))
    docker_args.append(str(threshold_i))
    
    # Run Docker container
    run_docker("docker.io/hedgelab/ciri_analysis", volumes, env_vars, docker_args)
    
    # Create results directory
    output_dir = os.path.join(main_mount_dir, "CIRI_assign_results")
    os.makedirs(output_dir, exist_ok=True)
    
    return Result(status="success", output_dir=output_dir)
  except Exception as e:
    logger.error(f"Docker execution failed: {str(e)}")
    return Result(status="error", output_dir="", message=str(e))


if __name__ == "__main__":
  import argparse
  
  parser = argparse.ArgumentParser(description="CRISPR QC and assignment analysis.")
  parser.add_argument('--target_directory', help="Working directory containing the h5 matrix and guides.csv")
  parser.add_argument('--script_directory', help="Directory containing the CIRI_unified.R script")
  parser.add_argument('--guide_strategy', choices=["1", "2"], help="Variable guides logic: 1=Single, 2=Dual")
  parser.add_argument('--matrix_filename', help="Filename of the HDF5 feature matrix")
  parser.add_argument('--threshold_a', help="Manual threshold for CRISPRa")
  parser.add_argument('--threshold_i', help="Manual threshold for CRISPRi")
  
  args = parser.parse_args()
  
  result = CIRI_assign(
    target_directory=args.target_directory,
    script_directory=args.script_directory,
    guide_strategy=args.guide_strategy,
    matrix_filename=args.matrix_filename,
    threshold_a=args.threshold_a,
    threshold_i=args.threshold_i,
  )
  
  print(f"Status: {result.status}")
  if result.status == "success":
    print(f"Output directory: {result.output_dir}")
  else:
    print(f"Error: {result.message}")
    sys.exit(1)
