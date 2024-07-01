#!/bin/sh
SDRranger \
count \
cDNA_fastqs \
--STAR-ref-dir=SDR001_REF_index \
--config=cDNA.json \
--output-dir=output_cDNA \
--threads=1 \
-vvv
