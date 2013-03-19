% Software Carpentry - Why Git?
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

## Provenance

Computational Scientist's motivation: reproduce, reproduce, reproduce!

## Reproduction

## Scientific Reproducibility

* **scientific provenance** - the origin or source of data presented
    in defense of a scientific theory or argument
* **experimental provenance** - equipment, materials, and procedure in
    obtaining data
* **computational provenance** - software, input data, and environment
    for obtaining data
**Never trust a result you can't reproduce**

## Managing Code Provenance

* **change management** How are changes allowed to be introduced?
* **release management** When do we (and what's needed to) publish
    this?
* **version management** Have we changed the way users interacted with
    our software?  Have we added functionality?  Changed it?  Removed
    it?

## Code Management

* **change control** Formal process for proposing, verifying, and
recording changes to software

### Computational Science Perspective

* not as relevant unless contributing to large or commercial software
engineering projects
* still must verify that scientific code produces the same results
after it has been changed

## Release Management

* **release management** - Management of the process of releasing
software to users or a production environment

### Computational Science Perspective:

* commercial/academic scientific software tools and libraries release
new versions of their code to users
* a researcher who publishes results is 'releasing' their code,
whether they realize it or not

## Release Management

software versioning schemes are somewhat arbitrary, generally:

* gcc-a.b.c

version numbers

* a (major) - functionality changed/removed from API
* b (minor) - functionality added to the API
* c (patch) - a bug was found and fixed


## Version Control
-------------------------------------

* **Version Control** - management of changes to documents, sources,
and other aspects a repository

### Computational Science Perspective **absolutely essential!**:

* aids and fosters collaboration
* encourages experimentation within repository
* removes needless commented out code cruft
* provenance for publishing scientific results


## History of Version Control
-------------------------------------

* 1972 - SCCS: commercial, introduced weave
* 1982 - RCS: open source, single file
* 1990 - CVS: open source, centralized repository
* 1998 - BitKeeper: commercial, distributed
* 2000 - Subversion: it's better than CVS!
* 2001 - Arch: open source, distributed
* 2003 - Monotone: open source, distributed, three-way merge
* 2005 - Hg: open source, fast, three-way merge
* 2005 - git: open source, very fast, many merge schemes

## GitHub
-------------------------------------

* offers free project hosting
   - so long as your code is public

* encourages collaboration and contributions

* similar to Google Code, Launchpad, BitBucket

Introduction to Git
-------------------------------------

* clone/init (start working on a repository)
* checkout - I made a mistake!
* status - what's different?
* diff - how is it different?
* add - files to a repository
* commit - changes to a repository

Looks very similar to SVN, but introduces no abstractions because
there is no central repository.
