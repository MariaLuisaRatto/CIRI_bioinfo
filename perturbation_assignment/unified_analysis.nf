// Nextflow Workflow: unified_analysis
// Unified CIRI workflow for Single and Dual guide analysis

// Input Parameters
params.scripts_dir = ''
params.data_dir = ''
// Allowed values: "1", "2"
params.mode = ''


process run_docker {
  container 'docker.io/hedgelab/ciri_analysis'
  input:
  val params.scripts_dir
  val params.data_dir
  val params.mode
  output:
  path 'results/'
  script:
    def args = [
    'Rscript',
    '/scripts/CIRI_unified.R',
    '/data',
    params.mode,
    ].join(' ')
    sh 'docker run --rm docker.io/hedgelab/ciri_analysis $args'
}
workflow {
  run_docker()
}
