# Documentation

[Back To
Testing](http://github.com/thehackerwithin/UofCSCBC2012/tree/master//) -
[Forward To Google Feedback
Form](https://docs.google.com/spreadsheet/viewform?formkey=dDlSWDEzMUt0Ri1TUDlTM21pUEwwSnc6MA#gid=0)

Just like version control and testing, documenting your code is the most
important thing you can do as a software developer. As we have seen in
previous sessions with other tools, good documentation is a sublime
experience that should permeate your code.

Documentation is important because it is [the only way that 90% of
people will ever interact with you or your
code](http://blip.tv/pycon-us-videos-2009-2010-2011/pycon-2011-writing-great-documentation-4899042).
In fact, it is the only way that scales up; there are only so many
emails that you can write.

What is disturbing is that documentation is a forgotten after thought
for most developers. It turns out that being able to write software and
being able to write in your primary spoken language are different
skills. Luckily, we are academics so it is in our nature to publish /
write. We have no excuse for bad documentation.

**The Many Stages of Documentation:**

1.  Readmes
2.  User Guides
3.  Developer Guides
4.  Code Comments
5.  API Documentation
6.  Auto-documentation

## Readmes

The omnipresent `README` file is typically a plaintext file that sits
next to the code. They typically may contain markup but are often quite
terse. The point of a readme file is to provide only the most basic of
information to the user / developer.

Readme files are so common that GitHub will render and display the
readme file for all directories whenever you are browsing a source tree.
Even Linux itself has a readme:

    Linux kernel release 3.x <http://kernel.org/>

> These are the release notes for Linux version 3. Read them carefully,
> as they tell you what this is all about, explain how to install the
> kernel, and what to do if something goes wrong.
>
> WHAT IS LINUX?
>
> > Linux is a clone of the operating system Unix, written from scratch
> > by Linus Torvalds with assistance from a loosely-knit team of
> > hackers across the Net. It aims towards POSIX and Single UNIX
> > Specification compliance.
> >
> > It has all the features you would expect in a modern fully-fledged
> > Unix, including true multitasking, virtual memory, shared libraries,
> > demand loading, shared copy-on-write executables, proper memory
> > management, and multistack networking including IPv4 and IPv6.
> >
> > It is distributed under the GNU General Public License - see the
> > accompanying COPYING file for more details.
>
> ...

## User's Guides

The next level of documentation are user's guides. These often take the
form of books or pdfs that aim to explain top level architechture and
functionality to possibly novice users. Such documents are extremely
helpful for bringing in new members to the community, going in depth
into the theory (math, biology, physics, chemistry, engineering), and as
a reference manual for advanced users and developers. However because of
their high level nature, you typically have to wait until the code has
stabalized to be able to write a good comprehensive user's guide.

**Examples:**
[FLASH](http://flash.uchicago.edu/site/flashcode/user_support/flash4b_ug.pdf),
[NumPy](http://www.tramy.us/numpybook.pdf).

## Developer Guides

Developer guides are very similar to user's guides except that they
assume a basic mastery of the project. They are typically for people who
want to *become* developers on a project rather than for existing
developers. They are probably most important for code projects that have
plugin architechtures and where the line between user and developer is
less well defined.

**Examples:** [Android](http://developer.android.com/guide/index.html),
[Python](http://docs.python.org/devguide/).
