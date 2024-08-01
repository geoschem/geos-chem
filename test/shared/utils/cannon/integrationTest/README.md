# Integration test utility scripts

This folder contains utility scripts that allow you to manually resubmit integration test jobs.  This can be useful, for example, if an integration test job fails but you don't want to re-create the entire directory structure or recompile the code.

## Scripts

`redoIntegrationTestCompile.sh`

  - Script to manually resubmit the integration test compilation job
  
  
`redoIntegrationTestExecute.sh`

  - Script to manually resubmit the integration test execution job
  
  
## Examples

1. Resubmit an integration test compilation job.

   ```console
   $ cd /path/to/integration/test/root/utils
   $ ./redoIntegrationTestCompile.sh
   ```

2. Resubmit an integration test execution job.

   ```console
   $ cd /path/to/integration/test/root/utils
   $ ./redoIntegrationTestExecute.sh
   ```

