#' Check if Docker is Available and Return Its Path
#'
#' The `has_docker` function checks if the Docker executable is available
#' in the system's `PATH`. It returns a list with a boolean indicating if Docker
#' is found, and the path to the Docker executable if it is available.
#'
#' @return
#' A list with two elements:
#' \item{found}{A logical value: `TRUE` if Docker is available in the system's
#' `PATH`, `FALSE` otherwise.}
#' \item{path}{A character string containing the full path to the Docker
#' executable if found,  or an empty string if not found.}
#'
#' @examples
#' result <- has_docker()
#' if (result$found) {
#'   cat("Docker is available at:", result$path, "\n")
#' } else {
#'   cat("Docker is not available.\n")
#' }
#'
#' @export
has_docker <- function() {
  path <- Sys.which("docker")
  return(list(found = nzchar(path), path = path, fn = run_in_docker))
}
#' Check if the script is running in Docker.
#'
#' @returns A truthy value indicating the state.
#' @export
is_running_in_docker <- function() {
  dockerenv_exists <- file.exists("/.dockerenv")
  cgroup_exists <- file.exists("/proc/1/cgroup")
  in_container_runtime <- FALSE
  if (cgroup_exists) {
    in_container_runtime <- any(
      grepl("docker", readLines("/proc/1/cgroup", warn = FALSE))
    )
  }
  return(dockerenv_exists || in_container_runtime)
}
#' Run a docker container.
#'
#' @param image_name The docker image you want to run.
#' @param volumes The list of volumes to mount to the container.
#' @param additional_arguments Vector of arguments to pass to the container.
#'
#' @export
run_in_docker <- function(image_name,
                          volumes = list(),
                          additional_arguments = c()) {
  base_command <- "run --privileged=true --platform linux/amd64 --rm"
  for (volume in volumes) {
    volume[1] <- normalizepath::normalize_path(volume[1],
      path_mappers = c(normalizepath::docker_mount_mapper)
    )
    base_command <- paste(base_command, "-v", paste(
      volume[1],
      volume[2],
      sep = ":"
    ))
  }
  base_command <- paste(base_command, image_name)
  for (argument in additional_arguments) {
    base_command <- paste(base_command, argument)
  }
  system2("docker", args = base_command, stdout = "", stderr = "")
}

#' CIRI_assign
#'
#' @description CRISPR QC and assignment analysis.
#' @param target_directory Working directory containing the h5 matrix and guides.csv
#' @param script_directory Directory containing the CIRI_unified.R script
#' @param guide_strategy Variable guides logic: 1=Single, 2=Dual (allowed values: 1, 2)
#' @param matrix_filename Filename of the HDF5 feature matrix
#' @param threshold_a Manual threshold for CRISPRa
#' @param threshold_i Manual threshold for CRISPRi
#' @return Results of the operation
#'
#' @export
CIRI_assign <- function(target_directory,
script_directory,
guide_strategy,
matrix_filename,
threshold_a,
threshold_i) {
  # Type validation
  if (!is.character(target_directory) || length(target_directory) != 1) {
    stop("target_directory must be a single character string")
  }
  if (!is.character(script_directory) || length(script_directory) != 1) {
    stop("script_directory must be a single character string")
  }
  valid_guide_strategy <- c("1", "2")
  if (!is.character(guide_strategy) || length(guide_strategy) != 1 || !(guide_strategy %in% valid_guide_strategy)) {
    stop(paste0("guide_strategy must be one of: ", paste(valid_guide_strategy, collapse=", ")))
  }
  if (!is.character(matrix_filename) || length(matrix_filename) != 1) {
    stop("matrix_filename must be a single character string")
  }
  if (!is.character(threshold_a) || length(threshold_a) != 1) {
    stop("threshold_a must be a single character string")
  }
  if (!is.character(threshold_i) || length(threshold_i) != 1) {
    stop("threshold_i must be a single character string")
  }
  
  # Security checks
  if (grepl("\\.\\./|\\.\\\\|\\/\\.\\./|\\\\\\.\\\\\\.\\\\", target_directory)) {
    stop("Path traversal detected in target_directory")
  }
  if (grepl("\\.\\./|\\.\\\\|\\/\\.\\./|\\\\\\.\\\\\\.\\\\", script_directory)) {
    stop("Path traversal detected in script_directory")
  }
  if (grepl("\\.\\./|\\.\\\\|\\/\\.\\./|\\\\\\.\\\\\\.\\\\", matrix_filename)) {
    stop("Path traversal detected in matrix_filename")
  }
  if (grepl("\\.\\./|\\.\\\\|\\/\\.\\./|\\\\\\.\\\\\\.\\\\", threshold_a)) {
    stop("Path traversal detected in threshold_a")
  }
  if (grepl("\\.\\./|\\.\\\\|\\/\\.\\./|\\\\\\.\\\\\\.\\\\", threshold_i)) {
    stop("Path traversal detected in threshold_i")
  }
  
  # Check if directory exists
  if (!is_running_in_docker()) {
    if (!dir.exists(target_directory)) {
      stop(paste("target_directory:", target_directory, "does not exist"))
    }
  }
  
  # Check if directory exists
  if (!is_running_in_docker()) {
    if (!dir.exists(script_directory)) {
      stop(paste("script_directory:", script_directory, "does not exist"))
    }
  }
  
  # Process file paths for Docker volume mounting
  # Process target_directory for Docker
  target_directory_abspath <- normalizePath(target_directory, mustWork = FALSE)
  target_directory_dir <- dirname(target_directory_abspath)
  target_directory_filename <- basename(target_directory)
  # Process script_directory for Docker
  script_directory_abspath <- normalizePath(script_directory, mustWork = FALSE)
  script_directory_dir <- dirname(script_directory_abspath)
  script_directory_filename <- basename(script_directory)
  
  # Main volume mount point
  main_mount_dir <- target_directory_dir
  
  # Execute Docker container with error handling
  tryCatch({
    result <- run_in_docker(
      image_name = "docker.io/hedgelab/ciri_analysis",
      volumes = list(
        c(target_directory_dir, "/mnt/analysis"),
        c(script_directory_dir, "/opt/scripts")
      ),
      additional_arguments = c(
        "/opt/scripts/CIRI_unified.R",
        "/mnt/analysis",
        guide_strategy,
        matrix_filename,
        threshold_a,
        threshold_i,
      )
    )
    
    # Process result
    return(list(
      status = "success",
      output_dir = file.path(main_mount_dir, "CIRI_assign_results")
    ))
  }, error = function(e) {
    stop(paste("Docker execution failed:", e$message))
  })
}
