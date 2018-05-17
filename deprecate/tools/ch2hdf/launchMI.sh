#!/bin/sh

LOCARGS="--block $COBALT_PARTNAME ${COBALT_CORNER:+--corner} $COBALT_CORNER ${COBALT_SHAPE:+--shape} $COBALT_SHAPE"
echo "Cobalt location args: $LOCARGS" >&2

   runjob $LOCARGS \
  --verbose=INFO \
  --label \
  --np 512 \
  --ranks-per-node 8 \
  --cwd $PWD \
  --exe $PWD/ch2hdf \
  --envs PAMI_DEVICE=B \
  --envs BG_MEMSIZE=16384 \
  --envs BG_THREADLAYOUT=1 \
  --envs XLSMPOPTS=parthds=1 \
  --envs PAMID_COLLECTIVES=1 \
  --envs PAMI_MEMORY_OPTIMIZED=1 \
  --envs BG_SHAREDMEMSIZE=512 \
  --envs BG_MAPCOMMONHEAP=0 \
  --args -simdata datawavelet000050.StreamerGridPointIterative.channel5 \
  --args -h5file output.channel5

#  --args -restart 0

#  --mapping TABCDE \
#  --envs BG_COREDUMPONEXIT=1 \

#  --envs BGLOCKLESSMPIO_F_TYPE=0x47504653


#  --envs PAMID_CONTEXT_POST=
#  --envs MUSPI_RECFIFOSIZE=67108864 \
#  --envs PAMID_CONTEXT_POST=1 \
