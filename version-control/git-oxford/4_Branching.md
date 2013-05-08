## A little on Branching

Up to this point we have been using what is called the **master** branch. Now, suppose you had heard of cascading style sheets and you wanted to do a little experimenting without affecting your web site. You can branch your development code to try something different - in a similar manner if you wanted to try some optimisation or add functionality without affecting your main version of your code you can create branches to try to do this. So lets try branching.

If you type:

     git branch
     * master

The "*" signifies what branch you are on and **master** is the only branch available to you. 

Now lets create a new branch:

     git branch css_test

Nothing appears to happen but when you type:

     git branch
       css_test
     * master

You will see that you have two branches and that you are on the **master** branch. Let's switch to the **css_test** branch. We do this by typing:

    git checkout css_test
    Switched to branch 'css_test'

If you do:

    git branch
    * css_test
      master

If you check what files you have in your repository you will see the same files as you had in your **master** branch. 

Let's now create a new css stylesheet, *mystyle.css* :

     /* make all paragraph text bold and red */
     p
     {
        color: red;
        font-weight: bold;
     }

Save this file, add it to git. Now link your style file
to your index.html file in the header section:

    <head>
      <title>My Home Page</title>
      <link rel="stylesheet" type="text/css" href="mystyle.css"/>
    </head>

If you now examine your index page on a browser you should see all the text appearing red and bold. Lets assume that you are happy with this state of affairs and you would like to add this effect to your web site. Add your index.html and commit both files.

Now go back to the master file:

    git checkout master
    Switched to branch 'master'

The first thing that you may notice is that your css file has now disappeared as have the changes to your html file - bummer! Rest assured they are still there but on **css_test** branch. You can go back to the **css_test** branch to check and come back to the **master** branch.

You are now happy that the world will want to see all your font bold and red so you want to incorporate it back into your **master** branch. You can do a merge:

    git merge css_test
    Fast-forward
    index.html  |    1 +
    mystyle.css |    7 +++++++
    2 files changed, 8 insertions(+), 0 deletions(-)
    create mode 100644 mystyle.css

Thankfully there were no conflicts. Now you will see the mystyle.css file and changes to your index.html in your master branch. It is useful to keep your **css_test** branch to do future testing with css elements. You can also merge changes from your master branch on to the **css_test**.

This is a very brief introduction to git branching but hopefully you can begin to appreciate the power. Branching in git is very easy to do and very fast.

Previous: [Collaborating with our colleagues](3_Collaboration.md) Next: [Conclusions and further information](5_Conclusion.md)