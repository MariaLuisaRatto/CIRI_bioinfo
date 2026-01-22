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

#' unified_analysis
#'
#' @description Unified CIRI workflow for Single and Dual guide analysis
#' @param scripts_dir Host directory containing the R scripts
#' @param data_dir Host directory containing the input data
#' @param mode Analysis Mode: 1=Single Guide, 2=Dual Guide (allowed values: 1, 2)
#' @return Results of the operation
#'
#' @export
unified_analysis <- function(scripts_dir,
data_dir,
mode) {
  # Type validation
  if (!is.character(scripts_dir) || length(scripts_dir) != 1) {
    stop("scripts_dir must be a single character string")
  }
  if (!is.character(data_dir) || length(data_dir) != 1) {
    stop("data_dir must be a single character string")
  }
  valid_mode <- c("1", "2")
  if (!is.character(mode) || length(mode) != 1 || !(mode %in% valid_mode)) {
    stop(paste0("mode must be one of: ", paste(valid_mode, collapse=", ")))
  }
  
  # Security checks
  if (grepl("\\.\\./|\\.\\\\|\\/\\.\\./|\\\\\\.\\\\\\.\\\\", scripts_dir)) {
    stop("Path traversal detected in scripts_dir")
  }
  if (grepl("\\.\\./|\\.\\\\|\\/\\.\\./|\\\\\\.\\\\\\.\\\\", data_dir)) {
    stop("Path traversal detected in data_dir")
  }
  
  # Check if directory exists
  if (!is_running_in_docker()) {
    if (!dir.exists(scripts_dir)) {
      stop(paste("scripts_dir:", scripts_dir, "does not exist"))
    }
  }
  
  # Check if directory exists
  if (!is_running_in_docker()) {
    if (!dir.exists(data_dir)) {
      stop(paste("data_dir:", data_dir, "does not exist"))
    }
  }
  
  # Process file paths for Docker volume mounting
  # Process scripts_dir for Docker
  scripts_dir_abspath <- normalizePath(scripts_dir, mustWork = FALSE)
  scripts_dir_dir <- dirname(scripts_dir_abspath)
  scripts_dir_filename <- basename(scripts_dir)
  # Process data_dir for Docker
  data_dir_abspath <- normalizePath(data_dir, mustWork = FALSE)
  data_dir_dir <- dirname(data_dir_abspath)
  data_dir_filename <- basename(data_dir)
  
  # Main volume mount point
  main_mount_dir <- scripts_dir_dir
  
  # Execute Docker container with error handling
  tryCatch({
    result <- run_in_docker(
      image_name = "docker.io/hedgelab/ciri_analysis",
      volumes = list(
        c(scripts_dir_dir, "/scripts"),
        c(data_dir_dir, "/data")
      ),
      additional_arguments = c(
        "Rscript",
        "/scripts/CIRI_unified.R",
        "/data",
        mode,
      )
    )
    
    # Process result
    return(list(
      status = "success",
      output_dir = file.path(main_mount_dir, "unified_analysis_results")
    ))
  }, error = function(e) {
    stop(paste("Docker execution failed:", e$message))
  })
}
