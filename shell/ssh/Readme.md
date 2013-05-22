# SSH

SSH is a protocol for securely sending data across an insecure network. An
important feature of ssh is the ability to perform public key cryptography.
We are about to take a little digression, but it is important because SSH is
built in to git at a very low level.

Public key cryptography consists of a public key, a private key, and an
algorithm called a cypher. Information can be combined with the public key
using the cypher in such a way that it appears like nonsense to anyone not
holding the private key. The only way to recover the originial information is
to use the cypher and the private key; using the cypher and the public key only
generates more nonsense-looking data. In this way, a person can use a public
key, encrypt data using the cypher + public key, and send the encrypted data
over the network without fear that someone will intercept the information.
In pseudo-code this looks something like:

    EncrpytedMessage = encrypt(PublicKey, Cypher, PlainTextMessage)
    PlainTextMessage = decrypt(PrivateKey, Cypher, EncryptedMessage)

Trying decode a message encrypted with the PublicKey using the PublicKey will
result in gibberish:

    Gibberish = decrypt(PublicKey, Cypher, EncryptedMessage)

This way, only a recipient holding the private key associated with that public
key can read the original message.

In addition to secure communication, public key cryptography can provide proof
of identity and data integrity: you can know that a message sent by someone
hasn't been altered. A person uses their private key and the cypher on some
data to create a hash, then sends the data and the resulting hash to someone
holding the public key. The person on the other end can take the public key +
hash and verify the data wasn't changed in transit.

## SSHing to another machine

SSH stands for Secure Shell, because it was originally designed in part as a
protocol for encrypted access to a remote machine.  This was meant as a
replacement for earlier programs such as rlogin that performed the same
function.  The idea is that when you "ssh to" a remote machine, it allows you
to remotely log into that machine and start a shell session, which is then
piped back to your local shell over an encrypted network connection.  You are,
in effect, communicating with a remote shell process wrapped inside your local
shell.

To ssh to a remote machine, simply enter:

    $ ssh example.com

This will then prompt you for your login password on the remote machine.

Note that because you have to log in to the remote machine to start a shell
session, you must have an account (and generally a home directory) on that
machine.  By default, ssh will assume that your username on the remote 
machine is the same as your local username.  If that is not the case you may
specify a different username using the following syntax:

    $ ssh username@example.com

If you are doing lots of work on remote machines, it becomes inconvenient to
type in your password every time.  However, if you have a private SSH key, and
the remote machine has a copy of your public key, we can use public key
cryptography to securely (and effortlessly on our part) identify ourselves to
the remote machine.

## Generating an SSH key pair

    $ cd ~/.ssh

It will likely say "no such file or directory."

    $ ssh-keygen -t rsa -C "your_email@youremail.com"
    Generating public/private rsa key pair.
    Enter file in which to save the key (/home/swc/.ssh/id_rsa):  <press enter>

The path that it provides will be to this home directory. This is okay. **Press
enter.** You may enter a passphrase. You'll see something like this :

    Created directory '/home/swc/.ssh'.
    Enter passphrase (empty for no passphrase):
    Enter same passphrase again:
    Your identification has been saved in /home/swc/.ssh/id_rsa.
    Your public key has been saved in /home/swc/.ssh/id_rsa.pub.
    The key fingerprint is:
    09:06:c6:0f:24:b7:84:ef:22:74:de:95:f0:99:64:5d your_email@youremail.com
    The key's randomart image is:
    +--[ RSA 2048]----+
    |  .+*   . .E     |
    |  .=o+ o .       |
    |   ..oB +        |
    | . ....B .       |
    |. o.. . S        |
    |. ....           |
    | . .             |
    |                 |
    |                 |
    +-----------------+

Now if we copy the content of our public key (id_rsa.pub) into a special file
called .ssh/authorized_keys in our home directory on the remote machine, we can
log in remotely without entering a password, so long as we have access to the
matching private key locally.

Having this set up will save us considerable hassle when we start wanting to
push to remote version control repositories later.


