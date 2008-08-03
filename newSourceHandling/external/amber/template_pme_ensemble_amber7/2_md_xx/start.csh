#!/usr/bin/env csh

set folder=`pwd`

echo "Starting $folder"

ssh %(hosts_prod)s "cd $folder; ./start_local.csh"
