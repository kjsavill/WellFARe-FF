#! /bin/bash

g09 < $1 > $1.out

grep "^ Energy= " $1.out | sed s/" Energy= "/""/ | sed s/"NIter=.*"/""/

