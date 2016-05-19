#!/bin/sh

dirname=$1

cd $1
ls | egrep '\.r8$' | xargs gzip -f
cd ..
