#!/bin/bash

NAME=freqMonoBit_AES.txt
rm $NAME
for (( i=1; i<=1024; i++ ))
do
  ./simul >> $NAME
done 
