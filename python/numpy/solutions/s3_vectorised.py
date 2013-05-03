import numpy as np
import operator

class Rec(object):
    def __init__(self, rawData):
        """Create a recommender from raw ratings.
        rawData is a dict mapping user names to their ratings.
        ratings are themselves dicts mapping paper ids to a rating.

        We create a (sparse) array self.ratings where the (i,j)th entry is the rating
        given by person i to paper j.
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

    def score(self, refIndex, ratings):
        """Score how similar each individual is to refIndex based on given ratings."""
        # Determine where individuals have rated the same items as refIndex
        mask = np.logical_and(ratings[refIndex,:] > 0, ratings > 0)
        # Calculate scores as the square root of the sum of squared rating differences
        return np.sqrt(((ratings[refIndex,:] - ratings)**2).sum(axis=1))

    def similarities(self, refIndex, ratings):
        """Find how similar each individual is to refIndex based on given ratings.

        This method takes ratings as an argument so that it can be used to rank
        both papers and people, just by transposing the ratings array.

        It returns a list of (index, score) pairs sorted in descending score order.
        """
        scores = self.score(refIndex, ratings)
        # Create & return the ordered list of pairs
        similarities = [(i, scores[i]) for i in range(len(scores)) if i != refIndex]
        similarities.sort(key=operator.itemgetter(1), reverse=True)
        return similarities

    def indexToName(self, pairs, names):
        """Convert indices in the first entry of each pair into names."""
        return map(lambda (idx, score): (names[idx], score), pairs)

    def similarPeople(self, person, n=3):
        """Find similar researchers based on ratings."""
        idx = self.people.index(person)
        similarities = self.indexToName(self.similarities(idx, self.ratings), self.people)
        return similarities[:n]
    
    def similarPapers(self, paper, n=3):
        """Find similar papers based on ratings."""
        idx = self.papers.index(paper)
        similarities = self.indexToName(self.similarities(idx, self.ratings.T), self.papers)
        return similarities[:n]

    def recommend(self, person, n=5):
        """Recommend papers for the given person to read, based on rating similarities."""
        idx = self.people.index(person)
        # Get the similarities scores for other researchers relative to this person
        scores = self.score(idx, self.ratings)
        # Create a mask over the entire ratings array that is only True where this person
        # hasn't rated a paper but the other person has
        mask = np.logical_and(self.ratings[idx, :] == 0, self.ratings > 0)
        # Compute the quantities upon which our recommendations are based
        totals = (self.ratings * mask * scores.reshape(7,1)).sum(axis=0)
        similarities = (mask * scores.reshape(7,1)).sum(axis=0)
        # Calculate a sorted list of recommendations
        recommendations = []
        for i, paper in enumerate(self.papers):
            if similarities[i] != 0:
                recommendations.append((paper, totals[i]/similarities[i]))
        recommendations.sort(key=operator.itemgetter(1), reverse=True)
        return recommendations[:n]
