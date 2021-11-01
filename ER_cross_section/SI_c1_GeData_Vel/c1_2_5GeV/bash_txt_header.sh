#!/bin/sh

awk  'BEGIN {print "double DM_2_5GeV_1000V[data_bin][2] = {";} { print "{ " $1",   " $2" }," } ENE { print "};"; }' DM_2_5GeV_10000V.txt     >  DM_2_5GeV_10000V.h 
awk  'BEGIN {print "double DM_2_5GeV_1167V[data_bin][2] = {";} { print "{ " $1",   " $2" }," } ENE { print "};"; }' DM_2_5GeV_11670V.txt     >	DM_2_5GeV_11670V.h 	
awk  'BEGIN {print "double DM_2_5GeV_1333V[data_bin][2] = {";} { print "{ " $1",   " $2" }," } ENE { print "};"; }' DM_2_5GeV_13330V.txt     >	DM_2_5GeV_13330V.h 	
awk  'BEGIN {print "double DM_2_5GeV_1500V[data_bin][2] = {";} { print "{ " $1",   " $2" }," } ENE { print "};"; }' DM_2_5GeV_15000V.txt     >	DM_2_5GeV_15000V.h 	
awk  'BEGIN {print "double DM_2_5GeV_1667V[data_bin][2] = {";} { print "{ " $1",   " $2" }," } ENE { print "};"; }' DM_2_5GeV_16670V.txt     >	DM_2_5GeV_16670V.h 	
awk  'BEGIN {print "double DM_2_5GeV_1883V[data_bin][2] = {";} { print "{ " $1",   " $2" }," } ENE { print "};"; }' DM_2_5GeV_18330V.txt     >	DM_2_5GeV_18330V.h 	
awk  'BEGIN {print "double DM_2_5GeV_2000V[data_bin][2] = {";} { print "{ " $1",   " $2" }," } ENE { print "};"; }' DM_2_5GeV_20000V.txt     >	DM_2_5GeV_20000V.h 	
awk  'BEGIN {print "double DM_2_5GeV_2333V[data_bin][2] = {";} { print "{ " $1",   " $2" }," } ENE { print "};"; }' DM_2_5GeV_23330V.txt     >	DM_2_5GeV_23330V.h 
awk  'BEGIN {print "double DM_2_5GeV_2167V[data_bin][2] = {";} { print "{ " $1",   " $2" }," } ENE { print "};"; }' DM_2_5GeV_21670V.txt     >  DM_2_5GeV_21670V.h 
awk  'BEGIN {print "double DM_2_5GeV_2500V[data_bin][2] = {";} { print "{ " $1",   " $2" }," } ENE { print "};"; }' DM_2_5GeV_25000V.txt     >  DM_2_5GeV_25000V.h 
awk  'BEGIN {print "double DM_2_5GeV_06667V[data_bin][2] = {";} { print "{ " $1",   " $2" }," } ENE { print "};"; }' DM_2_5GeV_06667V.txt    >	DM_2_5GeV_06667V.h
awk  'BEGIN {print "double DM_2_5GeV_08333V[data_bin][2] = {";} { print "{ " $1",   " $2" }," } ENE { print "};"; }' DM_2_5GeV_08333V.txt    >	DM_2_5GeV_08333V.h
awk  'BEGIN {print "double DM_2_5GeV_05000V[data_bin][2] = {";} { print "{ " $1",   " $2" }," } ENE { print "};"; }' DM_2_5GeV_05000V.txt    >	DM_2_5GeV_05000V.h
#awk  'BEGIN {print "double DM_2_5GeV_03333V[data_bin][2] = {";} { print "{ " $1",   " $2" }," } ENE { print "};"; }' DM_2_5GeV_03333V.txt    >	DM_2_5GeV_03333V.h







