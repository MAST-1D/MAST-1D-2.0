#!D:\pypy\pypy.exe
# EASY-INSTALL-ENTRY-SCRIPT: 'virtualenv==12.1.dev0','console_scripts','virtualenv-2.7'
__requires__ = 'virtualenv==12.1.dev0'
import sys
from pkg_resources import load_entry_point

if __name__ == '__main__':
    sys.exit(
        load_entry_point('virtualenv==12.1.dev0', 'console_scripts', 'virtualenv-2.7')()
    )
