Documentation
________________________________
`Back To Testing  <http://github.com/thehackerwithin/UofCSCBC2012/tree/master//>`_ - 
`Forward To Google Feedback Form <https://docs.google.com/spreadsheet/viewform?formkey=dDlSWDEzMUt0Ri1TUDlTM21pUEwwSnc6MA#gid=0>`_

Just like version control and testing, documenting your code is the most important thing
you can do as a software developer.  As we have seen in previous sessions with other tools, 
good documentation is a sublime experience that should permeate your code.

Documentation is important because it is `the only way that 90% of people will ever interact
with you or your code`_.  In fact, it is the only way that scales up; there are only so 
many emails that you can write.  

What is disturbing is that documentation is a forgotten after thought for most developers. 
It turns out that being able to write software and being able to write in your primary
spoken language are different skills.  Luckily, we are academics so it is in our 
nature to publish / write.  We have no excuse for bad documentation.

.. _the only way that 90% of people will ever interact with you or your code: http://blip.tv/pycon-us-videos-2009-2010-2011/pycon-2011-writing-great-documentation-4899042

**The Many Stages of Documentation:**

#. Readmes
#. User Guides
#. Developer Guides
#. Self-Documenting Code
#. Code Comments
#. API Documentation
#. Auto-documentation

Readmes
==========
The omnipresent ``README`` file is typically a plaintext file that sits next to
the code.  They typically may contain markup but are often quite terse.  The 
point of a readme file is to provide only the most basic of information to the 
user / developer.  

Readme files are so common that GitHub will render and display the readme file 
for all directories whenever you are browsing a source tree.  Even Linux itself
has a readme::

            Linux kernel release 3.x <http://kernel.org/>

    These are the release notes for Linux version 3.  Read them carefully,
    as they tell you what this is all about, explain how to install the
    kernel, and what to do if something goes wrong. 

    WHAT IS LINUX?

      Linux is a clone of the operating system Unix, written from scratch by
      Linus Torvalds with assistance from a loosely-knit team of hackers across
      the Net. It aims towards POSIX and Single UNIX Specification compliance.

      It has all the features you would expect in a modern fully-fledged Unix,
      including true multitasking, virtual memory, shared libraries, demand
      loading, shared copy-on-write executables, proper memory management,
      and multistack networking including IPv4 and IPv6.

      It is distributed under the GNU General Public License - see the
      accompanying COPYING file for more details. 

    ...


User's Guides
=============
The next level of documentation are user's guides.  These often take the form of 
books or pdfs that aim to explain top level architechture and functionality to 
possibly novice users.  Such documents are extremely helpful for bringing in new
members to the community, going in depth into the theory (math, biology, physics, 
chemistry, engineering), and as a reference manual for advanced users and 
developers.  However because of their high level nature, you typically have to wait 
until the code has stabalized to be able to write a good comprehensive user's guide.

**Examples:** `FLASH`_, `NumPy`_.

.. _FLASH: http://flash.uchicago.edu/site/flashcode/user_support/flash4b_ug.pdf
.. _NumPy: http://www.tramy.us/numpybook.pdf


Developer Guides
================
Developer guides are very similar to user's guides except that they assume a
basic mastery of the project.  They are typically for people who want to *become*
developers on a project rather than for existing developers.  They are probably 
most important for code projects that have plugin architechtures and where the
line between user and developer is less well defined.

**Examples:** `Android`_, `Python`_.

.. _Android: http://developer.android.com/guide/index.html
.. _Python: http://docs.python.org/devguide/


Self-Documenting Code
=====================
Much like in testing where you can simply write perfect code the first time, 
there is an analogous philosophy is documentation.  This is the philosophy of
`self-documenting code`_.  This ethos makes the claim that it is often 
possible to write code in such a way that new readers can understand what the
code does simply by reading it.  Therefore, no extra documentation is required.
It is all there in the code itself.

While there are obvious pitfalls with this approach (assumed knowledge on the 
reader's behalf, unavoidable complexities, etc) there are some merits.  By 
having meaningful naming conventions and structure it does become possible to
infer a lot about a code just by glancing at it.  Coding standards come from the
same desire to have readable software.  

However using this documentation strategy exclusively is *highly* inadvisable.

.. _self-documenting code: http://c2.com/cgi/wiki?SelfDocumentingCode


Code Comments
=============
Every language has a special character (or two) which indicate to the parser, 
compiler, or interpreter that whatever comes after or between these characters
should be ignored.  This allows the author to write annotate and explain the 
code that they are writing *right at the point that they are writing it!*  This 
is especially helpful if something wierd, obtuse, or obscure is about to happen
because it gives the author a chance to explain themselves to future developers
(often themselves in 1, 2, 6 months).

The best part is that you can put literally *anything* in comments: pubilcation 
citations, ASCII art, messages to lost loves, and threats to your collaborators.

In Python, the comment character is the hash symbol ``#``.  The following example
shows how you might help explain a toaster:

.. code-block:: python

    def toast(slices, toastiness, msg=None):
        # make sure the toaster has the right setting
        toastiness = int(toastiness) if 0 < toastiness else 5

        print "Engage the bread warming!"
        for slice if slices:
            slice.toast(toastiness)

        # log the message, making a default if needed
        if msg is None:
            msg = "Toasted to level {}".format(toastiness)
        logger.info(msg, "toast")

However, it is certainly possible to over-document your code with comments.  
Comments should never simply repeat what the code itself is doing.  The goal is to 
strike the right balance.  The approriate ratio changes with language.  (Typically
higher level languages have greater functionality per line and so have more comments.)
Try to avoid the following:

.. code-block:: python

    # init a to 0
    a = 0

    # make b 'a string'
    b = 'a string'

    # Add one to a
    a = a + 1

    # stopping excessive comments
    self.fall_on_sword()


API Documentation
=================
The application programming interface (API) is the definition of the protocol that 
two pieces of code may use to interact with one another.  Consider the case of 
functions.  All functions have a function signature which specifies how many 
arguments they accept and their return values.  This signature along with the 
module name and function name is the API.  (The function object/pointer itself 
is the implementation and is independent of the abstract API.)

Just because you have an argument list, however, does not imply that the meaning
of the arguments is known.  For example:

.. code-block:: python

    def f(a, b=10):
        ...

We know that ``f()`` accepts two argument ``a`` and ``b`` and that ``b`` should probably
be an integer.  But what does ``f()`` actually do?  What do these arguments mean in this
context?

Python allows the user to define API documentation right at the function, class, module, 
or variable definition.  Every Python object may have an ``__doc__`` attribute which is a string 
representation of the API docs.   This is known as a *docstring*.  `PEP257`_ describes
the convetions for docstrings.  The most important of these is that simple things should 
have simple docstrings.

Right below a definition, if the first non-comment, non-whitespace line is an 
unassigned string literal, then this string is automatically loaded in as the docstring.
It is this docstring which then read by the ``help()`` builtin or the ``?`` in IPython.

.. code-block:: python



.. _PEP257: http://www.python.org/dev/peps/pep-0257/
