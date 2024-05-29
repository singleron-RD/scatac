#!/usr/bin/env python

import sys

import snapatac2

bam, sample = sys.argv[1:3]
output_file = f"{sample}.fragment.gz"
snapatac2.pp.make_fragment_file(bam, output_file, barcode_regex="(\D*):\d*$")
