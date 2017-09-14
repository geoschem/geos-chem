#!/bin/bash

for a in $(ls *.f90); do mv -v $a ${a%.f90}.F90; done

exit 0
