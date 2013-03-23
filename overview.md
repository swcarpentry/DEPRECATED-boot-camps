% Software Carpentry - Overview
  American University of Beirut
% Aron Ahmadia
  <aron@ahmadia.net>
  <http://aron.ahmadia.net>
% 20 March, 2013

## Copy This Lecture!
<br></br>
<br></br>
<br></br>
<br></br>
<br></br>
<br></br>
<a rel="license" href="http://creativecommons.org/licenses/by/3.0/deed.en_US"><img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by/3.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" href="http://purl.org/dc/dcmitype/InteractiveResource" property="dct:title" rel="dct:type">Software Carpentry Overview</span> by <a xmlns:cc="http://creativecommons.org/ns#" href="http://aron.ahmadia.net" property="cc:attributionName" rel="cc:attributionURL">Aron Ahmadia</a> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/3.0/deed.en_US">Creative Commons Attribution 3.0 Unported License</a>.

## Introduction

Our goal:

<!-- There are a plethora of techniques, tools, languages, and
platforms available to help you scientifically compute. It is likely
that you will only be able to afford a limited amount of time learning
a subset of them. The purpose of this bootcamp is to help orient you
on the path to writing software as part of your research by:
-->

>
* introducing the important tools and paradigms of software construction with:
    + examples
    + hands-on exercies
* relating them to your role as a scientist, and
* making specific recommendations for
    + selecting tools that you should be using
    + and practices that you should be following

>
* three themes I will focus on are:
    + **reducing complexity**
    + **automating your workflow**
    + **encouraging experimentation**

Keywords in **bold** are important concepts, notes in *italics* are my own
personal advocacy.

# More About Software Carpentry

## History

* Founded by Greg Wilson in 1998.
* Open sourced materials 2004-06
* Currently funded by Sloan Foundation

## What We Teach

* Unix Command Line Interface (Shell)
* Version Control
* Python
* Testing

## What We *Actually* Teach

* A program is just another piece of lab equipment
* Programming is a human activity
* Little pieces loosely joined
* Let the computer repeat it
* Paranoia makes us productive
* Better algorithms beat better hardware

*How to THINK like a programmer*

## Who We Teach

<div align="center">
<table>
<tr>
<td><img src="swc-demographics/thumb-age.png" /></td>
<td><img src="swc-demographics/thumb-role.png" /></td>
</tr>
<tr>
<td><img src="swc-demographics/thumb-gender.png" /></td>
<td><img src="swc-demographics/thumb-platform.png" /></td>
</tr>
</table>
</div>

## My Goals for You

You will understand the principles of:

* Command Line Interface and the Bash Shell
* Python as a Scientific Computing Platform
* Collaborating with Git and GitHub
* Continuously Verifying and Validating Your Code

## Schedule

Today:

* Understanding the Shell
* Version Control with Git and Collaborating with GitHub
* Scientific Computing with Python (Start)

Tomorrow:

* Scientific Computing with Python (Finish)
* Verifying, Validating, and Driving Your Development with Tests
* Bring Your Code!

# Understand the Languages of Computing

## Be fluent in multiple languages
<comment>
You speak multiple languages when interacting with a computer.
Choosing to use a new tool, library, or language can be similar to
learning a new language:
</comment>

>
+ There is a high initial startup cost as you learn vocabulary, grammar, and
idioms
<font color=blue>`sum(x*y for x,y in itertools.izip(x_vector, y_vector))`
</font color>
+ You will learn faster by observing and working with others who are more
skilled than you
+ But once you have gained some fluency, you will find yourself capable of
new things!

## Use domain specific languages and libraries to increase your expressivity

* Aim for languages and tools that allow you to express your models simply,
whether that is symbolically, numerically, or
like [`chebfun`](http://www2.maths.ox.ac.uk/chebfun/),
a strange union of the two.

## Motivating Example: Chebfun

Chebfun teaser from the [website](http://www2.maths.ox.ac.uk/chebfun/)

```{.matlab}
% What's the integral of sin(sin(x)) from 0 to 10?
>> x = chebfun('x',[0 10]); sum(sin(sin(x)))

% What's the maximum of sin(x)+sin(x2) over the same interval?
>> x = chebfun('x',[0 10]); sum(sin(sin(x)))

% What's the solution to u"-xu=1 with zero boundary conditions on [-20,20]?
>> L = chebop(@(x,u)diff(u,2)-x.*u,[-20,20],'dirichlet'); plot(L\1)
```

## Motivating Example: FEniCS solving Poisson

[Dolfin Example 7](http://fenicsproject.org/documentation/dolfin/1.0.0/python/demo/pde/poisson/python/documentation.html), part of the [FEniCS project](http://fenicsproject.org/)

```{.python}
from dolfin import *
# Create mesh and define function space
mesh = UnitSquare(32, 32)
V = FunctionSpace(mesh, "Lagrange", 1)

# Define Dirichlet boundary (x = 0 or x = 1)
def boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

# Define boundary condition
u0 = Constant(0.0)
bc = DirichletBC(V, u0, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("10*exp(-(pow(x[0] - 0.5, 2) +
    pow(x[1] - 0.5, 2)) / 0.02)")
g = Expression("sin(5*x[0])")
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds

# Compute solution
u = Function(V)
solve(a == L, u, bc)
```

## Self-Promoting Example: PyClaw solving Acoustics

[PyClaw example](http://numerics.kaust.edu.sa/pyclaw/tutorial.html)
from [the paper submitted to the SIAM Journal of Scientific Computing](http://numerics.kaust.edu.sa/papers/pyclaw-sisc/pyclaw-sisc.html).

```{.python}

from clawpack import pyclaw
from clawpack import riemann

# Define the equations and the solution algorithm
solver = pyclaw.ClawSolver2D(riemann.rp2_euler_4wave)

# Define the problem domain
domain = pyclaw.Domain([0.,0.],[1.,1.],[100,100])
solver.all_bcs = pyclaw.BC.extrap
solution = pyclaw.Solution(solver.num_eqn,domain)

# Set physical parameters
gamma = 1.4
solution.problem_data['gamma']  = gamma

# Set initial data
xx,yy = domain.grid.p_centers
l = xx<0.5; r = xx>=0.5; b = yy<0.5; t = yy>=0.5
solution.q[0,...] = 2.*l*t + 1.*l*b + 1.*r*t + 3.*r*b
solution.q[1,...] = 0.75*t - 0.75*b
solution.q[2,...] = 0.5*l  - 0.5*r
solution.q[3,...] = 0.5*solution.q[0,...]* \
                    (solution.q[1,...]**2+solution.q[2,...]**2) \
                    + 1./(gamma-1.)

# Solve
solver.evolve_to_time(solution,tend=0.3)

# Plot
pyclaw.plot.interactive_plot()
```

## Use REPL Environments for Development

**REPL (read-eval-print-loop)** environments tighten the coupling between
the code you write and the results you see, increasing productivity.

REPL                                        non-REPL
-------------------          -----------------------
Visual Basic                 Supercomputers
Command Line                 C
IPython and Python           C++
MATLAB                       Java
Mathematica                  Fortran
LISP

*We will work with the command line and IPython for
the interactive lecture units today and tomorrow*

## Understand the limitations of floating-point numbers, learn symbolic tools as well
* As a computational scientist, you will likely be working with numerical
systems solved using floating point numbers
    + The excellent article *What Every Computer Scientist Should Know About
    Floating-Point Arithmetic* is a great way to do some advance reading

* Don't limit yourself to numerical computations, use symbolic tools such as
SAGE, SymPy and Mathematica as well when they make sense

# Don't Repeat Yourself (or Others)

## Automate common actions by saving simple blocks of code into **scripts**

>
* A script is a set of commands organized into a single file
* Sometimes it takes a few arguments, but more often there are just
some parameters at the top of the file to modify
* The script is the basest unit of scientific programming, you should be
comfortable writing these whenever you want to save or otherwise document or
repeat your actions
* Use scripts to explore new ideas, they are easy to write and throw away
* **Don't repeat commands into your REPL, save them to a script**

## Refactor commonly used blocks of code into **functions**
>
* Eventually, you will find that your scripts have a lot of repeated code,
or that you are spending a lot of time adjusting parameters at the top of
the file
* Refactor out repeated code into **function calls** in your scripts and
implement the **function** either in the same file or a separate one
* Be comfortable with the calling and return syntax of your programming language
environment, whether it is bash or Python
* **Don't repeat code in scripts, refactor them to functions**

## Group commonly used functions into **libraries**
>
* If you are unlucky enough to have to write a lot of software functions for
your work, you might want to consider designing and releasing a library so that
others do not have to share your misfortune
* You might want to first check that nobody else has implemented the
functionality you need
* If something close exists, it may be worth adapting to your needs if the
project is of high quality and suitably licensed
* *Openly licensed non-commercial libraries tend to have a much longer effective
lifespan than unreleased codes*
* **Share your code with others, and use their code**

## Design languages around commonly used paradigms

>
* I am not Guido van Rossum or Cleve Moler, but you could be!

# Reduce Complexity

## Basic strategies

>
* Endeavor to use languages and libraries that reduce the complexity of
your work
* It is worth installing a complicated or expensive software tool if your
computations are naturally expressed with it
* Always look for opportunities to write **less** code
    + you will have to do less initial work (sometimes)
    + you will introduce less bugs
    + your code will be easier to understand and maintain
* When writing software, try to keep individual functions short, single-purpose,
and avoid excessive nesting

# Organize, Checkpoint, and Collaborate

## Back up your data!

## Organize with wikis

* organize your personal notes into a personal wiki (gollum, gitit, instiki)
* organize your shared notes into a group wiki (gollum, gitit, instiki)
* *improve the layperson's understanding of your field by editing Wikipedia*

## Use version Control for checkpointing and collaboration

* use local version control software to checkpoint personal code development
  + checkpointing your work encourages wild ideas and late-night coding sessions
  + you can easily restore back in the morning if it was a bad idea
* use **distributed version control** to collaborate with others
* I advocate *git*, but you may be stuck with whatever your group uses
  + though check out git-svn for using git to collaborate with an svn repository,
  it's awesome!

*We will learn more about working with git and GitHub today and tomorrow*

# Verify and Validate your Code

## Principles of verification and validation
* **verification** - is your code correctly written?
* **validation** - do your computations accurately model the physical phenomena
in question?
* test frameworks help you verify your code, but validation is usually a manual
process
 + although it is desirable to write regression tests that verify previously
 validated results hold true when the code has been modified!
* use the **method of manufactured solutions** to verify correctness of code
* use **comparisons to experiment** to validate code
* use **comparisons to similar software** as an additional check

# Document your Computational Work

## Principles of documentation
* Save every bit of code you use for generating publishable results
* Document and comment your code for yourself as if you will need to understand
it in 6 months
  + use README files liberally
  + as well as file-level, function-level, and inline documentation
* If any piece of code is too complex to easily describe, consider refactoring
it

# Closing Thoughts

## Aim for reproducibility
* The goals of non-review scientific publications are to:
    + Describe a new result or approach
    + Convince others of its usefulness
* The **reproducibility** of your publication will greatly benefit both of
these goals
* See <http://
* See <http://figshare.com> for an easy way to store and share your data in
an easily-cited way

## Questions for the lecturer
* Why do senior computational scientists recommend against using libraries?
* How do I evaluate the usefulness of other people's code?

# References and Further Reading

# Books

## Pragmatic Programmer, The: From Journeyman to Master
Andrew Hunt and David Thomas

ISBN 978-0132119177

*Describes the important principles and practices of being an effective
programmer, instead of teaching a specific language or technique*

## Code Complete Second Edition
Steve McConnell

ISBN 978-0735619678

*Focuses on principles of software construction, with attention to skills,
testing, and design.*

## Verification and Validation in Scientific Computing
William L. Oberkampf and Christopher J. Roy

ISBN 978-0521113601

*Focuses on verification and validation of numerical solutions to models
described by systems of partial differential and integral equations.*

# Research Literature

## Programming Languages for Scientiﬁc Computing
Matthew G. Knepley

Preprint: http://arxiv.org/pdf/1209.1711.pdf

*Gives an overview of modern programming languages and techniques such as code
generation, templates, and mixed-language designs. This is a preprint,
so expect some rough spots.*

## Two Solitudes
Greg Wilson

Slides: http://www.slideshare.net/gvwilson/two-solitudes

*Describes Greg's journey as a scientist and leader for the Software Carpentry
project, provides some insight into the differences between industry and
academics.*

## Best Practices for Scientific Computing
D. A. Aruliah, C. Titus Brown, Neil P. Chue Hong, Matt Davis, Richard T. Guy,
Steven H. D. Haddock, Katy Huff, Ian Mitchell, Mark Plumbley, Ben Waugh,
Ethan P. White, Greg Wilson, Paul Wilson

Preprint: http://arxiv.org/abs/1210.0530

*Good summary paper of many fundamental practices for working with and
developing scientific software. This is a preprint, so expect some rough spots.*

# Web References

## What Every Computer Scientist Should Know About Floating-Point Arithmetic
David Golberg

Web article: http://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html

*Introduction to the IEEE floating-point standard, its implications, and many of
the common pitfalls when using floating-point numbers in scientific computing*

## The Research Software Engineer
Rob Baxter, Neil Chue Hong, Dirk Gorissen, James Hetherington, and Ilian Todorov

Web article: http://dirkgorissen.com/2012/09/13/the-research-software-engineer

*Discussion of the current challenges to scientific software engineering as
a profession.*

## Science Code Manifesto

http://sciencecodemanifesto.org

*Publicly signed commitment to clear licensing and curation of software
associated with research publications.*
