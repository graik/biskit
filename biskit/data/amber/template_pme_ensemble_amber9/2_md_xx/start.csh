#!/usr/bin/env csh

set folder=`pwd`

echo "Starting $folder"

ssh %(nodes_prod)s "cd $folder; ./start_local.csh"
