## Conclusions and further information

We've seen how we can use version control to,

* Keep track of changes like a lab notebook for code and documents.
* Roll back changes to any point in the history of changes to our files - "undo" and "redo" for files.
* Back up our entire history of changes in various locations.
* Work on our files from multiple locations.
* Identify and resolve conflicts when the same file is edited within two repositories without losing any work.
* Collaboratively work on code or documents or any other files.

Now, consider again our initial questions,

* Someone emails you, "Can I have a copy of the code you used to analyse the data that you graphed in a journal paper?". How easy would it be for you to get the code, or even the data, that generated that particular graph for them?
* If your laptop were stolen today how much would cost you to reproduce the work that went with it?
* How would you work on a piece of code, or a paper, on your home PC, a work PC and a laptop? How would you ensure that you kept all your files synchronized?
* How would you develop a piece of software collaboratively with colleagues? How would you ensure your sources were kept up to date? How would you deal with conflicts? 


Version control serves as a log book for your software and documents, ideas you've explored, fixes you've made, refactorings you've done, false paths you've explored - what was changed, who by, when and why - with a powerful undo and redo feature!

It also allows you to work with others on a project, whether that be writing code or papers, down to the level of individual files, without the risk of overwriting and losing each others work, and being able to record and understand who changed what, when, and why. 

> "If you are not using version control then, whatever else you may be doing with a computer, you are not doing science" -- Greg Wilson

### Find out more...

* K. Ram  (2013) "git can facilitate greater reproducibility and increased transparency in science", Source Code for Biology and Medicine 2013, 8:7 doi:[10.1186/1751-0473-8-7](http://dx.doi.org/10.1186/1751-0473-8-7) - survey of the range of ways in which version control can help research.
* [Visual Git Reference](http://marklodato.github.com/visual-git-guide/index-en.html) - pictorial representations of what Git commands do.
* [Pro Git](http://git-scm.com/book) - the "official" online Git book.
* [Version control by example](http://www.ericsink.com/vcbe/) - an acclaimed online book on version control by Eric Sink.
* [Git commit policies](http://osteele.com/posts/2008/05/commit-policies) - images on what Git commands to with reference to the working directory, staging area, local and remote repositories.
* [Gitolite](https://github.com/sitaramc/gitolite) - a way for you to host
your own multi-user Git repositories. Your collaborators send you their public SSH keys then they can pull and push from/to the repositories.
* G. Wilson, D. A. Aruliah, C. T. Brown, N. P. Chue Hong, M. Davis, R. T. Guy, S. H. D. Haddock, K. Huff, I. M. Mitchell, M. Plumbley, B. Waugh, E. P. White, P. Wilson (2012) "[Best Practices for Scientific Computing](http://arxiv.org/abs/1210.0530)", arXiv:1210.0530 [cs.MS].

Previous: [Collaborating with our colleagues](3_Collaboration.md) Next: [Git hints and tips](5_HintsAndTips.md)
