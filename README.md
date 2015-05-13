pb_calibration
==============

Filter and calibration programs for Illumina sequencing data (in BAM files).

(The name of this project no longer fits its functionality well so it may change...).

#Tools

## spatial_filter

Identify regions with spatially correlated errors e.g. bubbles, (from aligned BAM files where the read name can be parsed for spatial location) and allow filtering out or marking of the BAM fail bit for reads in those regions.

## calibration_pu & predictor_pb

Calculate and apply purity (calculated form dif files), cycle, and (good/bad) tile predictor based qualities using alignments and known SNPs file (typically calculated using alignments of spiked phiX and applied to the whole lane).


#Building

    cd src
    autoreconf --force --install
    ./configure --prefix=$PWD/.. --with-samtools=/software/solexa/pkg/samtools/current
    make install

