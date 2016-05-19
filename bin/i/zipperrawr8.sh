#!/bin/sh

dirname=$1

cd ./$dirname
ls | egrep '^im?p????*' | xargs gzip -f
cd ..
