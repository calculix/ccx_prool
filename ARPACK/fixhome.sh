#!/bin/sh
# this shell script is automatically fixed 'home' variable in file ARmake.inc
# prool ;-)
cp ARmake.inc ARmake.orig
echo "home="`pwd` >> ARmake.inc
