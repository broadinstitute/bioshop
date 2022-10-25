#!/usr/bin/env python

with open('JointVcfFiltering_input.wdl') as fh:
    wdl = fh.read()

indent = ' ' * 4
pyscript = ''

with open('python-script.py') as fh:
    for line in fh:
        line = indent + line
        pyscript += line

wdl = wdl % pyscript
with open('JointVcfFiltering_formatted.wdl', 'w') as fh:
    fh.write(wdl)
