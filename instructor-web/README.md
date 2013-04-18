# Instructor's notes for Finding Stuff Out From The Web

Google - obviously. How do we tell good from bad?
Stack Overflow
man pages / reference manuals / tutorials.
python documentation
pypi / cran / npm
github - also issues
mailing lists - archives
twitter

# Google

trial search:
 "github markdown to html"
 then I spotted that it's called "github flavoured markdown"
 and changed my search to
 "rendering github flavoured markdown"
 The top hit is the same. The second hit is more interesting.
 And the third hit is what I want.

Pasting bits of cryptic error messages into google often works
well. But when it doesn't work... xkcd.com/979

# Stack Overflow

It's a giant database of programming questions, and, more
importantly, answers. Anyone who is registered can ask a
question, answer a question, comment on an answer or a question,
and can upvote answers and questions. It seems to work pretty
well. 

Stack Overflow questions are often found near the top of a
google search, and when they are, I'll usually try those first.
Every question, answer, and comment on Stack Overflow is
referenceable as a URL and this is great for pasting into
comments to explain why the mysterious piece of code works the
way it does.

# Documentation

Almost every significant piece of software has documentation.
Some of it is even worth reading. I generally recommend that if
you're going to learn a language or a library you should follow
the tutorial if there is one, and at least skim the reference
documentation. Skimming gives you some idea of the landscape,
and hopefully when you're stuck later you'll remember "Ah,
didn't the Python Library Reference have something about CSV
files in it?".

For shell it's man bash. Which is huge, but again, I recommend
skimming, then dipping into on the beach, sipping your martini.

# Standards

Some languages, Fortran, C, JavaScript, POSIX shell, and software
technologies, SVG, PDF, are governed by standards. Which are dry technical
documents that tell you exactly what you can expect to work and
what may only work on certain implementations. Generally avoid
unless you have to, but they are the source of ultimate truth.

Be careful though, they are only standards. They may or may not
be followed to the letter. Don't do what I do and spend a day
implementing a feature that relied on the Trailing-Headers
feature of the HTTP specification, only to discover that no
web-browser or web-server implements Trailing-Headers. It's
perfectly well documented in the spec, it's just that it doesn't
work, and there's no notes telling me that.

# Ecosystems

Surrounding each programming language there is an ecosystem of
after-market modules. Usually there is a more or less standard
place to register modules and make them available. Python has
pypi, R has CRAN, perl has CPAN. These are places to find
modules, and also to find documentation.

# github

Source code for modules is found in various places, tarfiles on
ftp servers, sourceforge, googlecode, bitbucket, but
increasingly the common thing is to put the source code on
github. It has made coding social, and is the current
front-runner.

If you have the constitution for looking at other people's code
you can find all the answers you need. If you have enough time.
An increasing number of people are using github as the canonical
place to find documentation for a module. And it's social. If
you liked someone's image processing library, maybe you should
see what else they have written...

Every github project comes with an issue tracker, and you might
find that a ticket thas already been created for the issue
you're having. Maybe someone has suggested a workaround or
contributed a fix. Maybe you should.

# Mailing lists.

Communities form around pieces of software and that often means
a mailing list. Join the mailing list. Lurk. Be polite when
asking questions. If the mailing list has a public archive you
can search it, and you may find archived messages are sometimes
the result of a google search.

"numpy mailing list"

If the list is dauntingly large (the python mailing list), join
a smaller one! Maybe your instution has one?

Try and remember that everything you send to a list is basically
permanently recorded forever and searchable using google.
Naivety is acceptable, rudeness generally not.

# Web 2.0

The internet works because we all contribute to it. So if you've
found a useful answer on Stack Overflow and it works, give it
some love and upvote it! "how do i read a shape file in python".

If you've worked out how to do something tricky, and no-one else
seems to have documented, and you think it's a useful thing,
write a blog post about it.

# Twitter

Ask your mates. The modern version of which is to ask your
followers, on Twitter. I personally don't do this much, but it's
the done thing. You're unlikely to get definitive answers, but
you'll get leads.

# Time management.

A good general strategy is to allocate a budget to researching a
topic or question. The eXtreme Programming people call this
time-boxing. Spend a morning researching the best ridge
regression algorithm, or 2 hours trying to install TeX. At the
end of the allocated period use (or make a note of in your lab
book) the best solution you have found. If you've failed to find
anything, at least you bounded the amount of time spent, and
have avoided geting lost down a rabbit-hole. Switchings tasks
will improve your morale, and avoid you getting stuck in a
depressing loop trying to install the Python M2Crypto library
forever.
