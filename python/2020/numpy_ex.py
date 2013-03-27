import numpy as np
import math
import operator

class Rec(object):
    def __init__(self, rawData):
        """Create a recommender from raw ratings.
        rawData is a dict mapping user names to their ratings.
        ratings are themselves dicts mapping paper keys to a rating.
        """
        self.people = sorted(rawData.keys())
        papers = set()
        for ratings in rawData.itervalues():
            papers.update(ratings.keys())
        self.papers = sorted(papers)
        self.ratings = np.zeros((len(self.people), len(self.papers)), dtype=float)
        for person, ratings in rawData.iteritems():
            person_idx = self.people.index(person)
            for paper, rating in ratings.iteritems():
                paper_idx = self.papers.index(paper)
                self.ratings[person_idx, paper_idx] = rating
    
    def euclid(self, person1, person2, ord=2):
        """Compute the 2-norm between the ratings of the given individuals."""
        idx1 = self.people.index(person1)
        idx2 = self.people.index(person2)
        common_papers = (self.ratings[idx1] * self.ratings[idx2]) != 0 # Could use numpy.logical_and
        diffs = self.ratings[idx1,common_papers] - self.ratings[idx2,common_papers]
        #norm = math.sqrt((diffs**2).sum())
        return np.linalg.norm(diffs, ord=ord)
    
    def similarPeople(self, person, n=3):
        """Find similar researchers based on ratings."""
        similarity = [] # (person, score) pairs
        for other in self.people:
            if other != person:
                score = self.euclid(person, other, ord=None)
                similarity.append((other, score))
        similarity.sort(key=operator.itemgetter(1), reverse=True)
        #similarity.sort(key=lambda ps: ps[1], reverse=True)
        return similarity[:n]
    
    def similarPapers(self, paper, n=3):
        """Find similar papers based on ratings."""
        similarity = [] # (paper, score) pairs
        idx = self.papers.index(paper)
        for i, other in enumerate(self.papers):
            if other != paper:
                mask = np.logical_and(self.ratings[:,idx] != 0, self.ratings[:,i] != 0)
                score = np.linalg.norm(self.ratings[mask,idx] - self.ratings[mask,i])
                similarity.append((other, score))
        similarity.sort(key=operator.itemgetter(1), reverse=True)
        return similarity[:n]
