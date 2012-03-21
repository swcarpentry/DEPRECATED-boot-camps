#!/usr/bin/env python

# Viva!
# ~el Scopz

"""This utility converts rst to GitHub appropriate markdown.  Note that 
pandoc is not sufficient to have syntax highlighting work properly here.
"""
import os
import subprocess


def pandoc2gh(filename):
    with open(filename, 'r') as f:
        s = f.read()

    s = s.replace("~~~~ {.sourceCode .python}\n", "```python\n")
    s = s.replace("~~~~\n", "```\n")

    with open(filename, 'w') as f:
        f.write(s)


def main():
    # get file lists
    rstfiles = [os.path.join(p, f) for p, d, files in os.walk('.') for f in files]
    rstfiles = [f for f in rstfiles if f.endswith('.rst')]
    mdfiles = [f.rpartition('.')[0] + '.md' for f in rstfiles]

    # convert to md
    cmds = [['pandoc', rst, '-o', md] for rst, md in zip(rstfiles, mdfiles)]
    map(subprocess.check_call, cmds)

    # Munge pandoc md to gh md
    map(pandoc2gh, mdfiles)


if __name__ == '__main__':
    main()
