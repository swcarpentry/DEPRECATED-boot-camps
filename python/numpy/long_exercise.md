This exercise will use NumPy to create a recommendation engine for
academic papers.

The first step is to decide on recommendation criteria. It would seem
that we want to take into account how highly the paper was rated by
other people, but we also want to consider ratings that are from people
with similar backgrounds to the user of the system. This leads to
another criteria: the similarity between ratings across papers that both
users have already rated.

We will divide the code into three pieces.

First, we will take a list of previous ratings and store them in a NumPy
array.

Second, we will introduce measures of the similarity between two
papers or between two people's ratings.

Third, we will use those measures to generate recommendations for a
user.

For bonus points consider writing a class to solve this exercise, with
the various functions as methods of the class, and the data stored as
attributes.


Create a NumPy array with the recommendations
=============================================

Of all the millions of papers out there, most people have only read a
few. Since almost everyone has no opinion on almost every paper, the
data is very sparse. A good way to store sparse data is in a dictionary
for each user, where each ratings is stored as a unique paper identifier
and a rating.

We want to turn this data into an array. For this section, write a function
that takes the dictionary data as input and generates:

1.  A numpy array where element i,j is the rating of person i for paper
    j.

2.  A python list where element i is the name of person i.

3.  A python list where element j is the name of paper j.

Calculate Similarity
====================

We can think of each person's ratings as a real vector, so a measure of
the similarity of two researchers is just the distance between their
rating vectors. There are two ratings we will explore.

1. Add a function to compute the 2-norm (Euclidean distance) between two people's
   ratings. However, your function should only consider papers for which *both*
   people have provided a non-zero rating. If there are 4 papers, and Tommy has
   read papers 1,2, and 3 while Katy has read papers 1,3, and 4, then you should
   only compute the similarity using papers 1 and 3. If two people have no
   shared ratings then return 0.

2. (Optional) Add a function to compute the Pearson correlation between two
   vectors in the same way you computed the 2-norm above.

3. (Optional) Look up the Tanimoto distance function. Add a function to
   compute the Tanimoto distance.

Generate a Recommendation
=========================

There are a few ways to look at the recommendation data.

First, we could ask which researcher is most like you. Write a function
that takes a researcher id and identifies the 5 researchers whose
ratings are most like that researcher.

Second, we could ask which papers have the most similar ratings. Write a
function that takes a paper id and identifies the 5 papers whose ratings
are most like the paper. (Hint: could we reuse the code we've already
written and use the transpose function?)

Third, we could ask for recommended papers for a researcher. Write a
function that identifies the top 5 papers that a researcher should
read.  Keep in mind that the function should only return papers that
the researcher has *not* already rated. The recommendation score for
each paper could be calculated by summing the ratings given by other
researchers, weighted by the similarity between each researcher and
the subject. For each candidate paper j, compute this sum over
researchers i not equal to subj:

    sum(rating[i,j] * similarity[subj,i]) / sum(similarity[subj,i])

Look up Toby Segaran "Programming Collective Intelligence" for more
ideas.


Write some tests
================

Use Nose to write some tests for your code, extending the stub in
test_recommend.py. Specifically, think about how you can test the
input code and the distance measurements.

General Advice
==============

There are a few functions and concepts you will want to investigate in numpy.

1. numpy.logical\_and
2. numpy.cov
3. Broadcasting

Data
====

Some sample rating data is provided in the file data.py.  You can import
the data directly using:

```python
    from data import raw_scores
```
