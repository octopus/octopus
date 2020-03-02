#!/usr/bin/env bash
# 

cd `dirname "$0"`/..
if [ -x "$(which git)" ] && git log > /dev/null 2>&1 ; then
	git log --pretty=format:'%H' -n 1 | grep -v gpg
fi
