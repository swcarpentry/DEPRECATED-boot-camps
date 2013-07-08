"""Module to convert messages to and from Morse code"""

import string
import sys

class MorseTranslator:
    """This class can translate to and from morse code."""
    def __init__(self):
        self._letter_to_morse = {'a':'.-', 'b':'-...', 'c':'-.-.', 'd':'-..', 'e':'.', 'f':'..-.', 
                                 'g':'--.', 'h':'....', 'i':'..', 'j':'.---', 'k':'-.-', 'l':'.-..', 'm':'--', 
                                 'n':'-.', 'o':'---', 'p':'.--.', 'q':'--.-', 'r':'.-.', 's':'...', 't':'-',
                                 'u':'..-', 'v':'...-', 'w':'.--', 'x':'-..-', 'y':'-.--', 'z':'--..',
                                 '0':'-----', '1':'.----', '2':'..---', '3':'...--', '4':'....-',
                                 '5':'.....', '6':'-....', '7':'--...', '8':'---..', '9':'----.',
                                 ' ':'/' }

        self._morse_to_letter = {}

        for letter in self._letter_to_morse:
            morse = self._letter_to_morse[letter]
            self._morse_to_letter[morse] = letter

    def encode(self, message):
        """Encode the passed message into morse,
           and return the Morse code string"""
        morse = []

        for letter in message:
            letter = letter.lower()
            morse.append(self._letter_to_morse[letter])

        return string.join(morse," ")

    def decode(self, message):
        """Decode the passed Morse code message
           and return a string containing the decoded message"""

        english = []

        # Now we cannot read by letter. We know that morse letters are
        # separated by a space, so we split the morse string by spaces
        morse_letters = string.split(message, " ")

        for letter in morse_letters:
            english.append(self._morse_to_letter[letter])

        # Rejoin, but now we don't need to add any spaces
        return string.join(english,"")

if __name__ == "__main__":    

    translator = MorseTranslator()

    while True:
        print "Instruction (encode, decode, quit) :-> ",

        # Read a line from standard input
        line = sys.stdin.readline()
        line = line.rstrip()

        # the first line should be either "encode", "decode"
        # or "quit" to tell us what to do next...
        if line == "encode":
            # read the line to be encoded
            message = sys.stdin.readline().rstrip()

            print "Message is '%s'" % message
            print "Encoded is '%s'" % translator.encode(message)

        elif line == "decode":
            # read the morse to be decoded
            message = sys.stdin.readline().rstrip()

            print "Morse is   '%s'" % message
            print "Decoded is '%s'" % translator.decode(message)

        elif line == "quit":
            print "Exiting..."
            break

        else:
            print "Cannot understand '%s'. Instruction should be 'encode', 'decode' or 'quit'." % line

