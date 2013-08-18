#!/bin/bash

TYPES=( "real(DP)" "integer" "type(vector3d)")
TYPENames=( "Field_d" "Field_i" "Field_vec3d")
FileNameHints=( "Field_d" "Field_i" "Field_vec3d")
total=3
EXT="f90"
TEMPLATE_FILE="./list_mod.template.f90"

i=0
while test $i -lt $total
do
    FILEName="${FileNameHints[i]}_mod.${EXT}"
    echo "generating..${FILEName}"
    cpp -P -DListElemType="${TYPES[$i]}" -DListTypeName="${TYPENames[$i]}" $TEMPLATE_FILE $FILEName

  i=$((i+1))   
done