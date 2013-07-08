"""Module to encrypt and decrypt messages"""

import string
import sys

class Encryptor:
    """This class can encrypt and decrypt messages."""
    def __init__(self):
        self._letter_to_secret = {'a':'r', 'b':'x', 'c':'o', 'd':'j', 'e':'i', 'f':'t', 
                                  'g':'y', 'h':'b', 'i':'p', 'j':'3', 'k':' ', 'l':'u', 'm':'v', 
                                  'n':'s', 'o':'a', 'p':'q', 'q':'4', 'r':'h', 's':'k', 't':'w',
                                  'u':'0', 'v':'c', 'w':'7', 'x':'8', 'y':'5', 'z':'9',
                                  '0':'z', '1':'2', '2':'6', '3':'g', '4':'l',
                                  '5':'d', '6':'1', '7':'f', '8':'n', '9':'m',
                                  ' ':'e' }

        self._secret_to_letter = {}

        for letter in self._letter_to_secret:
            secret = self._letter_to_secret[letter]
            self._secret_to_letter[secret] = letter

    def encrypt(self, message):
        """Encrypt the passed message,
           and return the encrypted string"""
        encrypted = []

        for letter in message:
            letter = letter.lower()
            encrypted.append(self._letter_to_secret[letter])

        return string.join(encrypted,"")

    def decrypt(self, message):
        """Decrypt the passed encrypted code message
           and return a string containing the decrypted message"""

        english = []

        for letter in message:
            english.append(self._secret_to_letter[letter])

        return string.join(english,"")

if __name__ == "__main__":    

    encryptor = Encryptor()

    while True:
        print "Instruction (encrypt, decrypt, quit) :-> ",

        # Read a line from standard input
        line = sys.stdin.readline()
        line = line.rstrip()

        # the first line should be either "encrypt", "decrypt"
        # or "quit" to tell us what to do next...
        if line == "encrypt":
            # read the line to be encoded
            message = sys.stdin.readline().rstrip()

            print "Message is '%s'" % message
            print "Encrypted is '%s'" % encryptor.encrypt(message)

        elif line == "decrypt":
            # read the encoded message to be decoded
            message = sys.stdin.readline().rstrip()

            print "Encrypted is '%s'" % message
            print "Message is '%s'" % encryptor.decrypt(message)

        elif line == "quit":
            print "Exiting..."
            break

        else:
            print "Cannot understand '%s'. Instruction should be 'encrypt', 'decrypt' or 'quit'." % line
