#!/bin/bash
  echo "Running on GDB on node `hostname`"
xterm +bdc +cm -e gdb -x gdb_cmd --args $*
#gnome-terminal -e "gdb -x gdb_cmd --args "$*
exit 0
