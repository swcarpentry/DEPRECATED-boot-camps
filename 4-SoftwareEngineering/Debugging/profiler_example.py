@profile
def main():
    """A function that may take a while."""
    print "Welcome to the 80s!"
    print

    lyric1 = "You spin me right round,"
    lyric2 = "Right round like a record, baby!"

    # Spin some wheels
    n = 0
    while n < 6000:
        print lyric1
        n += 1

    # Spin some more wheels
    for n in range(4000):
        print lyric2

    print 
    print "Thanks to 'Dead or Alive' for that claaasic 1984 hit!"


# Run the function
if __name__ == "__main__":
    main()

