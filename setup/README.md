Students
========

This directory contains scripts for testing your machine to make sure
you have the software you'll need for your boot camp installed.  See
the comments at the head of each script for more details, but you'll
basically want to see something like:

    $ python swc-installation-test-1.py
    Passed
    $ python swc-installation-test-2.py
    check virtual-shell...  pass
    …
    Successes:

    virtual-shell Bourne Again Shell (bash) 4.2.37
    …

If you see something like:

    $ python swc-installation-test-2.py
    check virtual-shell...  fail
    …
    check for command line shell (virtual-shell) failed:
      command line shell (virtual-shell) requires at least one of the following dependencies
      For instructions on installing an up-to-date version, see
      http://software-carpentry.org/setup/
      causes:
      check for Bourne Again Shell (bash) failed:
        could not find 'bash' executable for Bourne Again Shell (bash)
        For instructions on installing an up-to-date version, see
        http://software-carpentry.org/setup/
    …

follow the suggestions to try and install any missing software.  For
additional troubleshooting information, you can use the `--verbose`
option:

    $ python swc-installation-test-2.py --verbose
    check virtual-shell...  fail
    …
    ==================
    System information
    ==================
    os.name            : posix
    …

