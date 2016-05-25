#!/bin/bash

DATATYPES=( "List_d" "List_vec3d" )
DATAELEMTYPES=( "real(DP)" "type(vector3d)" )
DATATYPE_HINTS=( "Scalar" "Vector" )
DATATYPEID=( 1 2 )
DATALOC=("vol" "surface" "point" )
EXT="f90"
TEMPLATE_FILE="./geometricfield_mod.template.f90"

for (( i=0; i < ${#DATATYPES[@]}; ++i ))
do
    for (( j=0; j < ${#DATALOC[@]}; ++j ))
    do
	FieldTypeName="${DATALOC[j]}${DATATYPE_HINTS[i]}Field"
	FILEName="${FieldTypeName}_mod.${EXT}"
	echo "generating..${FILEName}"
	cpp -P -DFieldTypeName="$FieldTypeName" -DFieldDataType="${DATATYPES[i]}" -DFieldDataElemType="${DATAELEMTYPES[i]}" -DDataDefLoc="DATADEFLOC_${DATALOC[j]}" -DDataTypeId="${DATATYPEID[i]}" $TEMPLATE_FILE $FILEName
	
    done
done
