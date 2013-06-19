UMass Amherst boot camp: Perspective from a helper
==================================================

Background
----------

A few weeks ago, I helped out with the [Software Carpentry boot camp at UMass Amherst](http://software-carpentry.org/bootcamps/2013-05-umass.html).  It was my first time being a helper.  Previously, my exposure to Software Carpentry was the MOOC style [online Spring 2011 class](http://software-carpentry.org/blog/2011/01/the-spring-2011-course-begins.html), which I took during physics graduate school.  I was among the 5-10% of students to finish the course.  The low completion rate was one of the reasons Software Carpentry switched to doing live, on-site boot camps.

A couple years later, I got in touch with Greg Wilson and he mentioned that there was a boot camp coming up in my area.  He invited me to join as a helper, I agreed, and I suddenly found myself on several mailing lists.  I had no idea what to expect, but the boot camp turned out to be quite fun and I met many interesting people -- including a core Python developer.

Post boot camp, I had many thoughts about the logistics and teaching.  The staff invited me to write these up in a blog post.

What went well
--------------

For the UMass Amherst boot camp, we had two instructors: Tracy Teal and Erik Bray.  They have very different backgrounds -- Tracy is a biological scientist, and Erik is a software engineer who works with physical scientists.  I felt this diversity was a real asset to the instruction because each instructor's teaching had a slightly different flavor.  Tracy emphasized the practical side of things more, whereas Erik would sometimes mention nuances that are intellectually interesting to a CS person.  Hopefully, that meant that there was something for everyone: practical stuff for the pragmatists, nuances for the computing nerds, variety for those who like variety. 

There was good energy in the room and people asked a lot of questions.  I liked how there were multiple ways to get help -- a) ask the teacher in front of the entire class, b) put a pink sticky on your laptop or raise your hand to get the attention of a helper, c) ask questions on the Etherpad. 

The [Etherpad](http://etherpad.org) was pretty cool, in my opinion, though I'm curious what the students thought.  Sometimes advanced questions were asked on the Etherpad, which was nice because the whole staff could see the questions and someone who understood that particular topic could chime in, and also because people could ask their expert or off-topic questions without intimidating beginners or disrupting the flow of the main instructor.  For example, someone wanted to know how to setup remote desktop software when the instructor was discussing SSH.  I even asked a couple questions on the Etherpad.

I experienced [IPython notebooks](http://ipython.org/notebook.html) for the first time, and that was neat.  I've used [Mathematica](http://www.wolfram.com/mathematica/) frequently for theoretical physics research, and it was nice to see a Mathematica notebook interface in Python.  My guess is that a student might have seen Python before, but they probably didn't know about IPython notebook.  So even if they already knew how to program in Python, at least the IPython notebook was novel. 

Technical difficulties
----------------------

I didn't expect to run into difficulties with software, since I was a helper and supposedly "good" at this.  Unfortunately, I didn't test all my software beforehand and that mistake haunted me.  The "cool", novel technologies I just mentioned gave me quite a few headaches.

Not everyone's IPython install worked, including mine.  I discovered that my version of IPython was 0.12 and not 0.13, since I had used apt-get to install IPython.  The Ubuntu IPython repository was not the latest.  I had also installed [Anaconda](http://continuum.io/downloads.html) just in case, but my path wasn't finding it, so I had to explicitly call the Anaconda IPython which was 0.13.  I was frantically trying to troubleshoot this problem during the first 15 minutes of Day 2.  I wish we had made everyone test their IPython install to make sure it was version >= 0.13, and that it properly loaded notebooks.  As far as I know, there is no explicit test for IPython right now in the Software Carpentry materials.

Even the Etherpad, which seems like a fairly foolproof thing (how can you mess up a web GUI tool?), was not without its problems.  As it turned out, my major responsibility at the boot camp was to take Etherpad notes on what the instructors were saying/doing.  I occasionally deleted some Etherpad text by accident and undo didn't restore everything.  I tried rewinding the history, but it wasn't ideal.  Later, I figured out that there is a save button and I tried to save more often.

It looks bad if the helpers can't even use the software properly.  I suggest in the future, to remind the helpers and instructors to test *all of their software*, especially tools they are unfamiliar with (e.g. IPython).  And tell whichever helper is taking Etherpad notes, to use the save button often!

Managing volunteers
-------------------

Overall, I think the staff could have done a better job managing the helpers and integrating them with the instructors.

Pre boot camp, I often had no idea what was going on.  The chatter on the mailing list was almost exclusively instructors discussing logistical matters with each other, which made me, as a helper, feel left out.  I had no idea how the waiting list worked (for the UMass boot camp, the number of spots was restricted) and I had no idea what I would be doing at the boot camp, until I explicitly asked about these things.  The organization could definitely have been improved.

During the boot camp, I didn't have as much interaction with the instructors and other helpers as I would have liked.  Any interaction that happened had to be initiated by myself -- for example, I hung around after the end of the sessions to chat.  Some of this was understandable since everyone was very busy and tired.  However, it would be nice if the instructors and helpers could hang out at some point before or after the boot camp -- have a meal, coffee, etc.

Managing volunteers is an important issue for Software Carpentry, in my opinion.  There are some organizations that take it very seriously, as evidenced by this [museum volunteer good practice checklist](https://docs.google.com/viewer?a=v&q=cache:IBDt8I5tWYkJ:www.mla.gov.uk/what/programmes/renaissance/regions/south_east/impact/~/media/South_East/Files/2010/V1%2520Museum%2520Volunteer%2520Management%2520Good%2520Practice%2520-%2520a%2520checklist.ashx+managing+museum+volunteers&hl=en&gl=us&pid=bl&srcid=ADGEESiEXMKLOwytyzOAEhAGNou8uoJxchCqK-AFyUr-v_nf0Gt67wf-QyEBcvce4qd45rjaS5FdHJss-QDg2ZMCNENz211sAzpabmh-uTPJpBnbmUce3rNTBE_yS9a-nNJNBjjFk7OA&sig=AHIEtbSmRh2WTQqXj76yr83R4Q0CaRO1EQ) I recently came across.

Software Carpentry only works if there are volunteers, so increasing the effort to make helpers feel welcome increases the probability that helpers will volunteer for future boot camps and spread good words about the organization.

Teaching Git
------------

I wasn't sure that people were understanding the material on Git.  Part of the problem is that Git is simply a tough subject to learn and teach, especially in such a limited time frame.  The lesson on Git started off with talking about SSH (because SSH keys are needed to interact with the GitHub repo) and then 5+ minutes of talking about abstract concepts in Git without any diagrams, slides, or visuals.  I was really confused about why we were spending so much time talking about SSH and not talking about Git.  I also think diagrams of the different repositories, commits, etc. could be introduced earlier.

On the plus side, the exercise of everyone pushing and pulling on the GitHub repo at the same time was pretty fun.  It definitely spiced up what could be a boring afternoon of version control.  More importantly, the exercise was both fun and educational.

Big picture and why scientists should care about software
---------------------------------------------------------

I've always thought that motivation is the most important part of giving a scientific talk, and the same should apply to boot camps.  That said, I'm not convinced we did the best job motivating our educational material at UMass.  I think it would be helpful to set aside a significant period of time to talk about why open and reproducible research is important, both at the beginning and the end of the boot camp.  I'm not sure that everyone even knows what open source software is, or the important role it's played in the history of mankind.  I can see that one wouldn't want to overly evangelize particular ideas and tools, but in my humble opinion, open source and reproducible research is worth a tiny bit of evangelizing.

In fact, this discussion could be fun and interesting.  One could show a few examples of scandals or big-time retractions and even ask the participants to share their experiences.  The recent headlines about the Excel error in austerity research comes to mind (testing anyone?)  The reason I got interested in Software Carpentry in the first place was because of a [Nature article](http://www.nature.com/news/2010/101013/full/467775a.html) about prominent scientific papers being retracted due to poor software.  I found that to be an attention-grabbing, compelling argument at the time.  (Greg was also quoted in the article.)  Show people both the carrot (you'll save so much time with software skills, people are more likely to cite you because your software works!) and the stick (career ending retraction of results, wasting your lab's time reproducing results that are wrong).

I know Software Carpentry tries to avoid PowerPoint style presentations, but I'm in favor of some slides on open and reproducible research.  As an example, I like the [wrap-up presentation](https://github.com/swcarpentry/boot-camps/tree/2013-06-southampton) that Robin Wilson wrote for the recent Southampton boot camp.  Robin reported that the [presentation went over pretty well](http://software-carpentry.org/blog/2013/06/soton-feedback.html).

Real-time assessment
--------------------

I wonder if anyone has tried in-class assessment at a boot camp, similar to the trend of using "clickers" in introductory science classes.  Angel Brady at Princeton wrote an excellent post [summarizing the ways in which you can do "virtual clickers"](http://blogs.princeton.edu/etc/2012/04/10/alternatives-to-physical-clickers-in-the-classroom).  I really have no idea how this would work for a Software Carpentry boot camp, but I'm just throwing the idea out there.  I am concerned that in-class assessment might be an extra layer of complication in a boot camp that is already stuffed with different technologies and material.

Following up after the boot camp
--------------------------------

Followup is important because it's so easy to be busy with other stuff after the boot camp is over and forget the material or lose energy and interest in using it.

### Office hours

Currently, Software Carpentry holds [office hours](http://software-carpentry.org/bootcamps/office-hours.html) after the boot camps.  This is a great idea, though I'm not sure how many people actually taken advantage of it, in practice.  Moreover, the office hours are pretty limited, no more than once a week.  Apparently, the last office hours were held in April, the UMass boot camp was held in May, and no new office hours have yet been announced on the webpage.

I think the office hours could be supplemented by other followup methods.

### Web tutorials

During the boot camp, I noticed a lot of interest on the Etherpad about supplemental material.  People were interested in tutorials for Python, Git, even vi/vim.  Instead of staff informally writing in some links on the Etherpad, maybe this could be done in a more organized fashion.  There might be a semi-official, curated list of internet resources.  Online software education has become a small industry in itself, and there are lots of great free web tutorials at places like [Code School](http://www.codeschool.com/paths/electives) and [Git Immersion](http://gitimmersion.com).

### IRC channels

One little known resource I want to emphasize is IRC channels.  They are certainly not mainstream; I've always been tech savvy and barely knew of their existence.  I thought they were a sketchy place to meet strangers.  It wasn't until I started doing open source software development that I used IRC and discovered what an outstanding resource it is.  I was able to ask questions on the #debian and #ubuntu channels and immediately receive incredibly helpful feedback on how to solve install problems.

I notice there are IRC channels for [R](http://www.r-bloggers.com/my-main-resources-for-r-programming/), [Python](http://www.python.org/community/irc/), [SciPy](http://www.scipy.org/scipylib/mailing-lists.html), and [Git](http://git-scm.com/community).  [Chatzilla](https://addons.mozilla.org/en-US/firefox/addon/chatzilla/) is a Firefox add-on and a nice IRC client that people could start with.  If people are intimidated about asking questions on IRC, someone more experienced with IRC can tag along as a buddy.  I would be happy to volunteer.

In fact, Software Carpentry has its own IRC channel #swcarpentry on [Freenode](http://freenode.net), and I suspect that most of the instructors don't even know about it.  At the moment, I'm working on open source software and am logged onto IRC all the time.  I'm happy to be available on #swcarpentry.  If Software Carpentry recruited more open source developers as helpers, they might be willing to staff the channel.

IRC channels have **huge potential** as a resource for Software Carpentry and I'd like to see them emphasized more at boot camps.  What could be better than a place online where you can interact with people real-time and get immediate help at almost any reasonable hour?

### Meetup groups

Further, if you live in the right area, there might be an amazing software related [Meetup](http://www.meetup.com) group near you.  [Bay Area R Users Group](http://www.meetup.com/R-Users/) and [Boston Python User Group](http://www.meetup.com/bostonpython/) come to mind.  I've personally attended a Bay Area R Users Group meetup and it was fantastic.  I met tons of interesting people, many of whom were scientists or had studied science in grad school.

### Emphasizing followup resources

For the Day 2 wrapup, there could be a brief "after SWC boot camp, where do you go from here?" speech outlining all the different and amazing resources available -- office hours, web tutorials, IRC, meetups, etc.  Even making people aware of meetup groups and IRC is an accomplishment in itself. 

Wrapping up
-----------

Volunteering at the boot camp brought my appreciation of Software Carpentry to a new level.  I already thought highly of the organization's work after taking the online class two years ago.  However, *actually being at a boot camp* impressed on me the enormous commitment and effort it takes to run a single boot camp -- then multiply this by dozens of boot camps that take place in a year.

I very much want to see Greg and the Software Carpentry team succeed.  I imagine a future where software skills are part of global science training, and where good software is universally accepted to be a **real and important part of scientific results**.  I will help as much as I can.
