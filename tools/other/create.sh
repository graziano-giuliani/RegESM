#!/bin/bash

ifile=$1

i1=`grep -n "INITINIT" $ifile | awk -F: '{print $1}'`
i2=`grep -n "ENDEND" $ifile | awk -F: '{print $1}'`
echo "       module mod_mit_gcm"
head -n $i2 $ifile | tail -n $((i2-i1))
echo "       contains" 
echo "       end module mod_mit_gcm"
