#!/bin/bash
set -eo pipefail

womtool validate JointVcfFiltering.wdl

fissfc meth_new \
    -d JointVcfFiltering.wdl \
    -s "Isolation forest joint VCF filtering" \
    -m JointVcfFiltering
