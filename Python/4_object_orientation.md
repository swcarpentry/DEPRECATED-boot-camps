
# Object Orientation

You have now learned how to package up your code into re-usable, documented functions, and how to then package up those functions into re-usable, documented modules (libraries). It is great to package up your code so that it is easy for other people to understand and re-use. However, one problem is that other people have a habit of re-using your code in the wrong, or in unexpected ways...

As an example, lets imagine someone using your [morse.py](3/example/morse.py) code using a module

    $ ipython
    $ import morse
    $ morse.letter_to_morse = "c"
    $ morse.encodeToMorse("Hello world")

    TypeError                                 Traceback (most recent call last)
    <ipython-input-5-8f98527ab404> in <module>()
    ----> 1 morse.encodeToMorse("Hello world")
    
         25     for letter in message:
         26         letter = letter.lower()
    ---> 27         morse.append(letter_to_morse[letter])
         28 
         29     return string.join(morse," ")

What happened here???

The problem is that the "letter_to_morse" variable is visible, and that we can change its value whenever we want. This breaks the "encodeToMorse" function, which expects "letter_to_morse" to be a dictionary.

As another example, lets imagine a user who is trying to edit the Morse code dictionary...

    $ ipython
    $ import morse
    $ morse.morse_to_letter["...."] = "K"
    $ morse.encodeToMorse("help")
    '.... . .-.. .--.'
    $ morse.decodeFromMorse(".... . .-.. .--.")
    'Kelp'

As you can see, allowing a user of your code to mess with the data on which it relies can lead to subtle, and difficult to find bugs!

## Encapsulation - Hiding Data

Object orientated programming solves this problem by packaging functions and their associated data together into a Class. A Class defines the type of data, together with functions that manipulate that data. Encapsulation is a key idea of object orientated programming, and means to hide the data in a Class, such that only the functions of defined as part of the Class can read or write (change) the data.

For example, take a look at this code from [guessgame.py](guessgame.py)

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

This piece of Python contains lots of new ideas. Before we explore them, lets try and play the game.

    $ ipython
    $ import guessgame
    $ game = guessgame.GuessGame("cat")
    $ game.guess("dog")
    Wrong answer. Try again!
    Out[3]: False
    $ game.guess("fish")
    Wrong answer. Try again!
    Out[3]: False
    $ game.guess("cat")
    Well done - you have won the game!
    Out[5]: True

Lets take a look at the help for GuessGame

    $ help(game)
    Help on instance of GuessGame in module guessgame:
    
    class GuessGame
     |  A simple guess the secret game
     |  
     |  Methods defined here:
     |  
     |  __init__(self, secret)
     |      Construct a game with the passed secret
     |  
     |  guess(self, value)
     |      See if the passed value is equal to the secret

"GuessGame", defined in this module is a example of a Class. Classes are used to package up functions with associated data. As you can see in the help(), we can only see the functions defined in the class. There are two functions, "__init__", which is used to construct a new Object of type GuessGame, and "guess" which is used to guess the secret. As you can see, the first argument to each of these functions is "self". "self" is a specialy variable that is used by the Class to gain access to the data hidden within.

Lets look again at the source for GuessGame (in [guessgame.py](guessgame.py))

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

Here you can see that the keyword "class" is used to define a new class (in this case, called GuessGame). Within the class you can see defined the two functions, "__init__" and "guess". The "__init__" function is special, and is called the "constructor". It must be present in all classes, and constructors are used in all object orientated programming languages. The job of the constructor is to define how to create an object of the class, i.e. how to initialise the data contained within an object instance of the class. In this case, you can see that the constructor specifies two variables, "_secret", which will hold the secret to be guessed, and "_nguesses", which holds the number of wrong guesses made to date. Note that these variables start with an underscore. This is the way you tell Python that the variables are private to the class. It is good programming practice to ensure that all class variable names in python are private, and start with an underscore.

Note that the variables are defined as attached to "self", via the full stop, e.g. "self._secret". "self" is a special variable that is only available within the functions of the class, and provides access to the hidden data of the class. You can see that "self" is used by the "guess" function to check the passed guess against "self._secret", and to increase the value of "self._nguesses" if the guess is wrong.

We can construct as many instances (objects) of a class as we want, and each will have its own "self" and its own set of hidden variables. For example;

    $ ipython
    $ from guessgame import GuessGame
    $ game1 = GuessGame("orange")
    $ game2 = GuessGame("carrot")
    $ game3 = GuessGame("apricot")
    $ game1.guess("apricot")
    Wrong answer. Try again!
    Out[4]: False
    $ game3.guess("apricot")
    Well done - you have won the game!
    Out[6]: True

(Note that we have used the "from X import Y" syntax in Python to import only GuessGame from [guessgame.py](guessgame.py). This allows us to write "game1 = GuessGame("orange")" instead of "game1 = guessgame.GuessGame("orange")".)

Note that we don't need to pass "self" ourselves to the class functions. "self" is passed implicitly by Python when we construct an object of the class, or when we call a function of the object.

## Exercise

### Exercise 4

Edit your [morse.py](3/example/morse.py) script create a class "MorseTranslator" by packaging together the functions "encodeToMorse" and "decodeFromMorse" with the variables "letter_to_morse" and "morse_to_letter".

Make sure that you document your class, e.g. by documenting the "__init__" function you will have to write, and also by documenting the class, as in the above GuessGame class in [guessgame.py](guessgame.py).

When you have finished, test that the Morse code produced by your class is correctly translated back to English, e.g.

    $ ipython
    $ from morse import MorseTranslator
    $ translator = MorseTranslator()
    $ message = "hello world"
    $ translator.decode( translator.encode(message) ) == message
    True

If you get really stuck, you can take a look at the completed example in [4/example/morse.py](4/example/morse.py).

### Extension

The act of encoding and decoding a message to and from Morse code is very similar to encrypting and decrypting a message. Try to write a new class, "Encryptor", that can encrypt a message using an "encrypt" function, and decrypt the message using a "decrypt" function. 

When you have finished, test that your Encryptor can decrypt its own encrypted messages, e.g.

    $ ipython
    $ from encryptor import Encryptor
    $ encryptor = Encryptor()
    $ message = "hello world"
    $ encryptor.decrypt( encryptor.encrypt(message) ) == message
    True

If you get really stuck, then you can take a look at the completed example in [4/example/encryptor.py](4/example/encryptor.py)

Make sure that you commit your edited script to your Git repository.

    $ git commit -am "...commit message..."
    $ git push

# [Previous](3_documenting_code.md) [Up](python_and_good_programming_practice.md) [Next](5_matplotlib.md) 
