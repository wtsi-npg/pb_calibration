pb_calibration
==============

Filter and calibration programs for Illumina BAM files

Building:
I've been unable to build this on anything other than lenny-dev64

cd src
./bootstrap
./configure --with-samtools=/software/solexa/bin/aligners/samtools/current
make

