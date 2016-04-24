#!/bin/bash

argstring=`./build_argstring.sh $@`

matlab -nodesktop -nosplash -r "try; standalone_peakdetection($argstring); catch ME; disp(getReport(ME,'extended','hyperlinks','off')); exit; end; exit"
