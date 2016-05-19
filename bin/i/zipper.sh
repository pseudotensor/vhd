#!/bin/sh

dirname=$1

cd $1
ls | egrep '\.ppm$' | xargs gzip
cd ..
