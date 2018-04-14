#!/bin/bash
set -ex
cd "$(bazel info workspace)"

bazel build cmd/bio-bam-sort/... cmd/bio-pamtool/...

rm -f /scratch-nvme/sortshard/*.sortshard

INBAM=/scratch-nvme/bam/CNVS-NORM-110033752-cfDNA-WGBS-Rep1.bam
OUTBAM=/scratch-nvme/tmp/foo.bam
OUTPAM=/scratch-nvme/tmp/foo.pam

bazel-bin/cmd/bio-bam-sort/sorter/go_default_test --test.run xxxx -test.v -test.bench Split -logtostderr -bam $INBAM --split
exit 1

bazel-bin/cmd/bio-bam-sort/bio-bam-sort -profile-interval-s 60 -heap-profile heap -cpu-profile cpu --bam $OUTBAM /scratch-nvme/sortshard/*.sortshard
sambamba index $OUTBAM

bazel-bin/cmd/bio-bam-sort/bio-bam-sort -profile-interval-s 60 -heap-profile heap -cpu-profile cpu --pam $OUTPAM /scratch-nvme/sortshard/*.sortshard

bazel-bin/cmd/bio-pamtool/bio-pamtool checksum $INBAM >/tmp/inbam.csum
bazel-bin/cmd/bio-pamtool/bio-pamtool checksum $OUTBAM >/tmp/outbam.csum
bazel-bin/cmd/bio-pamtool/bio-pamtool checksum $OUTPAM >/tmp/outpam.csum
