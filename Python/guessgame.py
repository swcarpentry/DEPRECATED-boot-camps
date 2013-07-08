"""A simple guessing game"""

class GuessGame:
    """A simple guess the secret game"""
    def __init__(self, secret):
        """Construct a game with the passed secret"""
        self._secret = secret
        self._nguesses = 0   # the number of guesses

    def guess(self, value):
        """See if the passed value is equal to the secret"""

        if self._nguesses >= 3:
            print( "Sorry, you have run out of guesses." )

        elif value == self._secret:
            print( "Well done - you have won the game!" )
            return True
        else:
            print( "Wrong answer. Try again!" )
            self._nguesses += 1  # increase the number of wrong guesses
            return False
