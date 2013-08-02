
letter_to_morse = {'a':'.-', 'b':'-...', 'c':'-.-.', 'd':'-..', 'e':'.', 'f':'..-.', 
                   'g':'--.', 'h':'....', 'i':'..', 'j':'.---', 'k':'-.-', 'l':'.-..', 'm':'--', 
                   'n':'-.', 'o':'---', 'p':'.--.', 'q':'--.-', 'r':'.-.', 's':'...', 't':'-',
                   'u':'..-', 'v':'...-', 'w':'.--', 'x':'-..-', 'y':'-.--', 'z':'--..',
                   '0':'-----', '1':'.----', '2':'..---', '3':'...--', '4':'....-',
                   '5':'.....', '6':'-....', '7':'--...', '8':'---..', '9':'----.',
                   ' ':'/' }

# We need to invert the dictionary. This will create a dictionary
# that can go from the morse back to the letter
morse_to_letter = {}

for letter in letter_to_morse:
    morse = letter_to_morse[letter]
    morse_to_letter[morse] = letter

message = "... --- ... / .-- . / .... .- ...- . / .... .. - / .- -. / .. -.-. . -... . .-. --. / .- -. -.. / -. . . -.. / .... . .-.. .--. / --.- ..- .. -.-. -.- .-.. -.--"

english = []

# Now we cannot read by letter. We know that morse letters are
# separated by a space, so we split the morse string by spaces
morse_letters = message.split(" ")

for letter in morse_letters:
    english.append(morse_to_letter[letter])

# Rejoin, but now we don't need to add any spaces
import string
print( string.join(english,"") )
