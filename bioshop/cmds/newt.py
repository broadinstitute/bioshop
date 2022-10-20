#!/usr/bin/env python

import argparse

from .. version import __VERSION__
from . import etl 
from . import fit 
from . import call 

Commands = ['etl', 'fit', 'call']


BANNER_TEXT = f"""
 _   _         _
| |_|_|___ ___| |_ ___ ___
| . | | . |_ -|   | . | . |
|___|_|___|___|_|_|___|  _|  
 {__VERSION__[:20]:<21}|_|
"""

def banner():
    msg = str.join('\n', filter(bool, BANNER_TEXT.split('\n')))
    msg += '\n'
    print(msg)

parser = argparse.ArgumentParser(prog='newt', description='Newt command runner')
subs = parser.add_subparsers(dest='command')
etl.get_cli_parser(subs)
fit.get_cli_parser(subs)
call.get_cli_parser(subs)

def main_cli():
    banner()
    args = parser.parse_args()
    if args.command == 'etl':
        return etl.main_cli(args=args)
    if args.command == 'fit':
        return fit.main_cli(args=args)
    if args.command == 'call':
        return call.main_cli(args=args)

if __name__ == "__main__":
    main_cli()
