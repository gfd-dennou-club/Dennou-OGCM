#!/bin/bash

ElemTYPES=( "real(DP)" "real(DP)" )
ElemTYPESIGS=( "d" "d" )
DIMS=( "2" "3" )
total=2
EXT="f90"
VECSPACE_TEMPLATE_FILE="./vectorspace.template.f90"
DOTPRODOPTR_TEMPLATE_FILE="./dotProd.template.f90"

i=0
while test $i -lt $total
do
    VecSpaceTypeName="vector${DIMS[$i]}${ElemTYPESIGS[$i]}"
    FILEName="${VecSpaceTypeName}_mod.${EXT}"
    echo "generating vector space .. ${FILEName}"
    cpp -P -Dvecspace_elem_type="${ElemTYPES[$i]}" -Dvecspace_elem_type_sig="${ElemTYPESIGS[$i]}" -Dvecspace_elem_size="${DIMS[$i]}" $VECSPACE_TEMPLATE_FILE $FILEName
    
    echo "generating dot product operator.. "
    dotProdFILEName="dotProd_${VecSpaceTypeName}_mod.${EXT}"
    cpp -P -DvectorspaceTypeName="${VecSpaceTypeName}" $DOTPRODOPTR_TEMPLATE_FILE $dotProdFILEName

    i=$((i+1))   
done