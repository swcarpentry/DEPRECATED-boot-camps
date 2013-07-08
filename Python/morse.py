
import string
import sys

letter_to_morse = {'a':'.-', 'b':'-...', 'c':'-.-.', 'd':'-..', 'e':'.', 'f':'..-.', 
                   'g':'--.', 'h':'....', 'i':'..', 'j':'.---', 'k':'-.-', 'l':'.-..', 'm':'--', 
                   'n':'-.', 'o':'---', 'p':'.--.', 'q':'--.-', 'r':'.-.', 's':'...', 't':'-',
                   'u':'..-', 'v':'...-', 'w':'.--', 'x':'-..-', 'y':'-.--', 'z':'--..',
                   '0':'-----', '1':'.----', '2':'..---', '3':'...--', '4':'....-',
                   '5':'.....', '6':'-....', '7':'--...', '8':'---..', '9':'----.',
                   ' ':'/' }

morse_to_letter = {}

for letter in letter_to_morse:
    morse = letter_to_morse[letter]
    morse_to_letter[morse] = letter


while True:
    print( "Instruction (encode, decode, quit) :-> ", )

    # Read a line from standard input
    line = sys.stdin.readline()
    line = line.rstrip()

    # the first line should be either "encode", "decode"
    # or "quit" to tell us what to do next...
    if line == "encode":
        # read the line to be encoded
        message = sys.stdin.readline().rstrip()

        print( "Message is '%s'" % message )
        print( "Encoded is '%s'" % encodeToMorse(message) )

    elif line == "decode":
        # read the morse to be decoded
        message = sys.stdin.readline().rstrip()

        print( "Morse is   '%s'" % message )
        print( "Decoded is '%s'" % decodeFromMorse(message) )

    elif line == "quit":
        print( "Exiting...")
        break

    else:
        print( "Cannot understand '%s'. Instruction should be 'encode', 'decode' or 'quit'." % line )

