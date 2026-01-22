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

def unified_analysis(scripts_dir: str, data_dir: str, mode: str) -> Result:
  """
  Unified CIRI workflow for Single and Dual guide analysis
  
  Parameters:
      scripts_dir: Host directory containing the R scripts
      data_dir: Host directory containing the input data
      mode: Analysis Mode: 1=Single Guide, 2=Dual Guide (allowed values: 1, 2)
  
  Returns:
      Result: Results of the operation
  """
  # Parameter validation
  if not isinstance(scripts_dir, str):
    raise TypeError(f"scripts_dir must be a string, got {type(scripts_dir).__name__}")
  scripts_dir_path = validate_path(scripts_dir)
  if not isinstance(data_dir, str):
    raise TypeError(f"data_dir must be a string, got {type(data_dir).__name__}")
  data_dir_path = validate_path(data_dir)
  mode_valid_values = ["1", "2"]
  if not isinstance(mode, str):
    raise TypeError(f"mode must be a string, got {type(mode).__name__}")
  if mode not in mode_valid_values:
    raise ValueError(f"mode must be one of {mode_valid_values}")
  
  # Directory existence checks
  if not is_running_in_docker():
    if not os.path.isdir(scripts_dir_path):
      raise NotADirectoryError(f"Directory {scripts_dir_path} does not exist")
  if not is_running_in_docker():
    if not os.path.isdir(data_dir_path):
      raise NotADirectoryError(f"Directory {data_dir_path} does not exist")
  
  # Process file paths for Docker volume mounting
  # Process scripts_dir for Docker
  scripts_dir_abspath = os.path.abspath(scripts_dir_path if 'scripts_dir_path' in locals() else scripts_dir)
  scripts_dir_dir = os.path.dirname(scripts_dir_abspath)
  scripts_dir_filename = os.path.basename(scripts_dir)
  # Process data_dir for Docker
  data_dir_abspath = os.path.abspath(data_dir_path if 'data_dir_path' in locals() else data_dir)
  data_dir_dir = os.path.dirname(data_dir_abspath)
  data_dir_filename = os.path.basename(data_dir)
  
  # Main volume mount point
  main_mount_dir = scripts_dir_dir
  
  # Execute Docker container with error handling
  try:
    # Prepare Docker volumes
    volumes = {}
    volumes[scripts_dir_dir] = "/scripts"
    volumes[data_dir_dir] = "/data"
    
    # Prepare environment variables
    env_vars = {}
    
    # Prepare Docker arguments
    docker_args = []
    docker_args.append("Rscript")
    docker_args.append("/scripts/CIRI_unified.R")
    docker_args.append("/data")
    docker_args.append(str(mode))
    
    # Run Docker container
    run_docker("docker.io/hedgelab/ciri_analysis", volumes, env_vars, docker_args)
    
    # Create results directory
    output_dir = os.path.join(main_mount_dir, "unified_analysis_results")
    os.makedirs(output_dir, exist_ok=True)
    
    return Result(status="success", output_dir=output_dir)
  except Exception as e:
    logger.error(f"Docker execution failed: {str(e)}")
    return Result(status="error", output_dir="", message=str(e))


if __name__ == "__main__":
  import argparse
  
  parser = argparse.ArgumentParser(description="Unified CIRI workflow for Single and Dual guide analysis")
  parser.add_argument('--scripts_dir', help="Host directory containing the R scripts")
  parser.add_argument('--data_dir', help="Host directory containing the input data")
  parser.add_argument('--mode', choices=["1", "2"], help="Analysis Mode: 1=Single Guide, 2=Dual Guide")
  
  args = parser.parse_args()
  
  result = unified_analysis(
    scripts_dir=args.scripts_dir,
    data_dir=args.data_dir,
    mode=args.mode,
  )
  
  print(f"Status: {result.status}")
  if result.status == "success":
    print(f"Output directory: {result.output_dir}")
  else:
    print(f"Error: {result.message}")
    sys.exit(1)
