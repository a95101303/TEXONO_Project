#!/bin/sh

sed  -i 's/D/E/g' *.txt 
sed  -i 's/\}/\n/g'  *.txt
sed  -i 's/{/ /g'  *.txt
sed  -i 's/,/ /g'  *.txt

