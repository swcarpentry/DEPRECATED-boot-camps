
letter_to_morse = {'a':'.-', 'b':'-...', 'c':'-.-.', 'd':'-..', 'e':'.', 'f':'..-.', 
                   'g':'--.', 'h':'....', 'i':'..', 'j':'.---', 'k':'-.-', 'l':'.-..', 'm':'--', 
                   'n':'-.', 'o':'---', 'p':'.--.', 'q':'--.-', 'r':'.-.', 's':'...', 't':'-',
                   'u':'..-', 'v':'...-', 'w':'.--', 'x':'-..-', 'y':'-.--', 'z':'--..',
                   '0':'-----', '1':'.----', '2':'..---', '3':'...--', '4':'....-',
                   '5':'.....', '6':'-....', '7':'--...', '8':'---..', '9':'----.',
                   ' ':'/' }

message = "SOS We have hit an iceberg and need help quickly"

morse = []

for letter in message:
    letter = letter.lower()
    morse.append(letter_to_morse[letter])

# We have the letters in a list. The "string" module can join this
# list back into a string, using "string.join". The first argument
# is the list, the second is the separater. In this case, we will
# join the list with spaces

import string
print( string.join(morse," ") )

