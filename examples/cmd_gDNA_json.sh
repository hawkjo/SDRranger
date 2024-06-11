#!/bin/sh
SDRranger \
count_gDNA \
gDNA_fastqs \
--STAR-ref-dir=SDR001_REF_index \
--config=gDNA.json.gz \
--output-dir=output_gDNA \
--threads=1 \
-vvv
