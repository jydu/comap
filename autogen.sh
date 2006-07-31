#! /bin/sh
touch NEWS README AUTHORS ChangeLog
aclocal
autoconf
automake --add-missing --copy
