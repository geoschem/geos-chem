# Parallel test utility scripts

This folder contains utility scripts that allow you to manually resubmit parallel test jobs.  This can be useful, for example, if an parallel test job fails but you don't want to re-create the entire directory structure or recompile the code.

## Scripts

`redoParallelTestCompile.sh`

  - Script to manually resubmit the parallel test compilation job
  
  
`redoParallelTestExecute.sh`

  - Script to manually resubmit the parallel test execution job
  
  
## Examples

1. Resubmit an parallel test compilation job.

   ```console
   $ cd /path/to/parallel/test/root/utils
   $ ./redoParallelTestCompile.sh
   ```

2. Resubmit an parallel test execution job.

   ```console
   $ cd /path/to/parallel/test/root/utils
   $ ./redoParallelTestExecute.sh
   ```

