## Resubmission
This is a resubmission for V1.0.0 of scDiffCom.

I have made the following changes:
* I have added data.table::setDTthreads(2) in the tests to avoid the NOTEs about parallel processing

## Test environments
* local ubuntu 22.04.2 install, R 4.3.1
* windows-latest, release, github action
* macOS-latest, release, github action
* ubuntu-latest, release, github action
* ubuntu-latest, devel, github action
* ubuntu-latest, oldrel-1, github action
* windows, devel, https://win-builder.r-project.org/

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs

## Downstream dependencies
There are currently no downstream dependencies for this package
