#!/bin/bash
argstring=""
i=0
for var
do
    if [ $i == 0 ]
    then
        if echo "$var" | egrep -q '^\-?[0-9]*\.?[0-9]+$'
        then
            argstring=$var
        else
            argstring=\'$var\'
        fi
    else
        if echo "$var" | egrep -q '^\-?[0-9]*\.?[0-9]+$'
        then
            argstring=$argstring,$var
        else
            argstring=$argstring,\'$var\'
        fi
    fi
    ((i++))
done
echo "$argstring"
