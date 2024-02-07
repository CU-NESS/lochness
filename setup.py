#!/usr/bin/env python
"""
File: setup.py
Author: Joshua J. Hibbard
Date: Feb 2024.

Description: Installs LOCHNESS.
"""
import os
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name='lochness', version='0.1',\
    description='Lunar Observatory Code in Healpy by the NESS team',\
    author='Joshua J. Hibbard',\
    author_email='joshua.hibbard@colorado.edu',\
    packages=['lochness'])

LOCHNESS_env = os.getenv('LOCHNESS')
cwd = os.getcwd()
if not LOCHNESS_env:
    import re
    shell = os.getenv('SHELL')
    print("\n")
    print("#" * 78)
    print("It would be in your best interest to set an environment variable")
    print("pointing to this directory.\n")
    if shell:
        if re.search('bash', shell):
            print("Looks like you're using bash, so add the following to " +\
                "your .bashrc:")
            print("\n    export LOCHNESS={0}".format(cwd))
        elif re.search('csh', shell):
            print("Looks like you're using csh, so add the following to " +\
                "your .cshrc:")
            print("\n    setenv LOCHNESS {!s}".format(cwd))
    print("\nGood luck!")
    print("#" * 78)
    print("\n")
elif LOCHNESS_env != cwd:
    print("\n")
    print("#" * 78)
    print("It looks like you've already got an lochness environment variable " +\
        "set but it's \npointing to a different directory:")
    print("\n    LOCHNESS={!s}".format(lochness_env))
    print("\nHowever, we're currently in {!s}.\n".format(cwd))
    print("Is this a different lochness install (might not cause problems), or " +\
        "perhaps just")
    print("a typo in your environment variable?")
    print("#" * 78)
    print("\n")
