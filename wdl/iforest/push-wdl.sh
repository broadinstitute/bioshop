#!/bin/bash
set -eo pipefail

./format_wdl.py 
womtool validate JointVcfFiltering_formatted.wdl

fissfc meth_new \
    -d JointVcfFiltering_formatted.wdl \
    -s "Isolation forest joint VCF filtering" \
    -m JointVcfFiltering
