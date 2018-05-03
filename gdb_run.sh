#!/bin/bash
echo "Running on GDB on node `hostname`"
export GTEST_CATCH_EXCEPTIONS=0
xterm +bdc +cm -e gdb -x cmd_gdb --args $*
exit 0
