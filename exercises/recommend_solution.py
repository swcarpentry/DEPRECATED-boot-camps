import numpy


EPS = 1.0e-9


def prep_data():
    """Reads the data file ratings.txt and returns an NxM array where each
    row is a person, and each column is that person's rating of a paper.
    Each person and paper is given a numeric ID associated with their sort
    order.
    """

    data = numpy.genfromtxt('../data/ratings.txt', dtype=None, names=True,
                            delimiter=', ')
    people = sorted(numpy.unique(data['Name']))
    papers = sorted(numpy.unique(data['Paper']))
    ratings = numpy.zeros((len(people), len(papers)))

    # Fill the ratings array
    for person, paper, rating in data:
        # Note: This would be very ineffcient over a large data set.  Why?
        # What sort of data structure would you use instead?
        person_id = people.index(person)
        paper_id = papers.index(paper)
        ratings[person_id, paper_id] = rating

    return people, papers, ratings


def top_matches(ratings, person, num, similarity_func):
    """Return the top n most similar individuals to a given person."""

    # Note that the similarity_func is actually a function that we passed in as
    # an argument--in this case either similarity_distance or
    # similarity_pearson.  We can do this in Python, and because the functions
    # take the same arguments and return the same *type* of output we can
    # expect this to work.

    scores = []
    for other in xrange(ratings.shape[0]):
        if other != person:
            scores.append((similarity_func(ratings, person, other), other))

    scores.sort()
    scores.reverse()
    return scores[:num]


def calculate_similar(paper_ids, ratings, num=10):
    """
    Find the papers that are most similar to each other.
    """

    result = {}

    # Now we want to look at the ratings by paper--that is, each row contains
    # the ratings for a given paper.
    ratings_by_paper = ratings.T

    for item in xrange(ratings_by_paper.shape[0]):
        # Use similarity_disance to find the top-scoring matches for a given
        # paper
        unnamed_scores = top_matches(ratings_by_paper, item, num,
                                     similarity_distance)

        scores = [(x[0], paper_ids[x[1]]) for x in unnamed_scores]
        result[paper_ids[item]] = scores

    return result


def similarity_distance(prefs, left_index, right_index):
    """Computes the similarity of two sets of ratings using the sum of squares
    metric:

        similarity(A, B) = 1 / (1 + norm(A, B))

    where

        norm(A, B) = distance**2 = (A1 - B1)**2 + ... + (An - Bn)**2
    """

    # Where do both people have preferences? Recall that the prefs array
    # contains one row per person, where the row data contains their ratings
    # of each paper (one column per paper).
    #
    # Recall that each person has not rated every paper sample.  Create a
    # boolean *mask* representing the papers each person has rated.  That is,
    # a 1 x N array where N is the number of papers, where each item is True
    # if the person has rated that paper, and False if they have not.  In
    # this data set 0 represents no rating (not a zero rating).

    # First, fill in the '???' to select the rows from prefs for each person's
    # ratings:
    left_prefs = prefs[left_index, :]
    right_prefs = prefs[right_index, :]

    # Then create the mask for which papers have been rated by each person:
    left_has_prefs = left_prefs > 0
    right_has_prefs = right_prefs > 0

    # Now combine the two masks into a single mask of papers rated by *both*
    # people using a logical AND (there's a function in the numpy module for
    # this that we saw earlier):
    mask = numpy.logical_and(left_has_prefs, right_has_prefs)

    # Return zero if there are no common ratings:
    if mask.sum() == 0:
        return 0

    # Now let's compute the sum of squares. First we need to get the delta
    # of each element in the vector arrays, while applying the mask to select
    # only the non-zero ratings
    diff = left_prefs[mask] - right_prefs[mask]

    # Use numpy.linalg.norm to compute the actual difference, but remember we
    # have to *square* it to get the sum of squares:
    sum_of_squares = numpy.linalg.norm(diff) ** 2

    # Recall, the final formula for the similarities is 1 / (1 + norm(A, B)):
    return 1.0 / (1.0 + sum_of_squares)


def similarity_pearson(prefs, left_index, right_index):
    """Computes the similarity of two sets of ratings using Pearson's
    correlation metric:

        similarity(A, B) = cov(A, B) / (stdev(A) * stdev(B))

    using numpy.cov() to compute the covariance matrix:

       var(X)  | cov(X, Y)
     ---------------------
     cov(X, Y) |  var(Y)

    where stdev(X) == sqrt(var(X))
    """

    # Create the mask where both people have ratings.  This is exactly
    # the same as in similarity_distance().
    # BONUS: Try squashing this down to as few lines as possible (try 3)
    left_prefs = prefs[left_index, :]
    right_prefs = prefs[right_index, :]

    mask = numpy.logical_and(left_prefs > 0, right_prefs > 0)

    # Note that summing an array of booleans gives the number of trues since
    # True == 1, False == 0
    if numpy.sum(mask) == 0:
        return 0

    # Calculate the Pearson score r
    # Use the Numpy function numpy.cov() mentioned in the docstring to compute
    # the var-covar matrix for the two rating sets
    varcovar = numpy.cov(left_prefs[mask], right_prefs[mask])

    # Extract the covariance from the matrix in the format shown in the
    # docstring.  This is the numerator in the formula for r
    numerator = varcovar[0, 1]

    # Extract the variances var(A) and var(B) from the matrix and use
    # numpy.sqrt() to compute the standard deviations for both rating sets
    denominator = numpy.sqrt(varcovar[0, 0]) * numpy.sqrt(varcovar[1, 1])

    # Why would we want to be extra sure the denominator is not close to zero?
    if denominator < EPS:
        return 0

    return numerator / denominator


def recommend(prefs, subject, similarity_func):
    """
    Get recommendations for an invdividual from a weighted average of other
    people.
    """

    totals = {}
    sim_sums = {}

    # Using the prefs array (remember rows are people) use the .shape attribute
    # of Numpy arrays to get the number of people and the number of papers in
    # the array
    num_people = prefs.shape[0]
    num_papers = prefs.shape[1]

    for other in xrange(num_people):

        # Don't compare people to themselves
        if other == subject:
            continue

        # Call the similarity function with the correct arguments.  Look back
        # at similarity_distance and similarity_pearson to check their function
        # signatures if you don't remember.
        similarity = similarity_func(prefs, subject, other)

        # Ignore scores of zero or lower
        if similarity < EPS:
            continue

        for title in xrange(num_papers):
            # Only score papers this perosn hasn't seen yet
            if prefs[subject, title] > EPS or prefs[other, title] < EPS:
                continue

            # Similarity * score
            if title not in totals:
                # Set the *default* similarity score for this title to 0
                totals[title] = 0
            else:
                # Add the similarity score for this title to the existing total
                # for that tile
                # Look up the rating by other(person), paper title and
                # multiply it by the similarity score to add to the total
                totals[title] += prefs[other, title] * similarity

            # Sum of the similarity scores
            if title not in sim_sums:
                sim_sums[title] = 0
            else:
                sim_sums += similarity

    # Create the normalized list

    rankings = []
    for title, total in totals.items():
        rankings.append((total/sim_sums[title], title))

    # Return the sorted list
    rankings.sort()
    rankings.reverse()
    return rankings


def test():
    person_ids, paper_ids, all_ratings = prep_data()
    print 'person_ids', person_ids
    print 'paper_ids', paper_ids
    print 'all_ratings', all_ratings
    print 'similarity distance', similarity_distance(all_ratings, 0, 1)
    print 'similarity Pearson', similarity_pearson(all_ratings, 0, 1)
    print top_matches(all_ratings, 0, 5, similarity_pearson)
    print calculate_similar(paper_ids, all_ratings)
    print recommend(all_ratings, 0, similarity_distance)
    print recommend(all_ratings, 1, similarity_distance)


if __name__ == '__main__':
    test()
