#!/bin/bash -x

SCALE=$1
cd /root/eriedy3/spgemm/files
unxz --keep spgemm-${SCALE}.bin.xz
emu_multinode_exec 3600 -- ../GrB-mxm-timer-baseline --verbose=100 -k 1,2,4 -f spgemm-${SCALE}.bin --binary -s ${SCALE}
rm spgemm-${SCALE}.bin

