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
    left_prefs = prefs[???]
    right_prefs = prefs[???]

    # Then create the mask for which papers have been rated by each person:
    left_has_prefs = left_prefs ???
    right_has_prefs = right_prefs ???

    # Now combine the two masks into a single mask of papers rated by *both*
    # people using a logical AND (there's a function in the numpy module for
    # this that we saw earlier):
    mask = numpy.???

    # Return zero if there are no common ratings:
    if mask.sum() == 0:
        return 0

    # Now let's compute the sum of squares. First we need to get the delta
    # of each element in the vector arrays, while applying the mask to select
    # only the non-zero ratings
    diff = left_prefs[???] ???

    # Use numpy.linalg.norm to compute the actual difference, but remember we
    # have to *square* it to get the sum of squares:
    sum_of_squares = numpy.linalg.norm(???) ???
    
    # Recall, the final formula for the similarities is 1 / (1 + norm(A, B)):
    return ???


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
    left_prefs = ???
    right_prefs = ???

    left_has_prefs = ???
    right_has_prefs = ???

    mask = ???
