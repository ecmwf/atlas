#!/usr/bin/env bash
# runtests.sh --
#    Bourne shell script to control a program that uses funit
#    Name of the program: first argument
#
#    $Id: runtests.sh,v 1.2 2008/01/26 11:15:10 arjenmarkus Exp $
#
if test -f ftnunit.lst ; then
    rm ftnunit.lst
fi
if test -f runtests.log ; then
    rm runtests.log
fi
echo ALL >ftnunit.run

RETVAL=0
chk=1
until test ! -f ftnunit.lst -a $chk -eq 0 ; do
    chk=0
    $(readlink -e $1) $2 $3 $4 $5 $6 $7 $8 $9
    RETVAL=$?
done

rm ftnunit.run

exit $RETVAL