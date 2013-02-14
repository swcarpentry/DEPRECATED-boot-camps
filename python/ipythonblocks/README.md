# Practicing Python with `ipythonblocks`

## Learning Goals

- **Part 1:** Color
    - Use IPython's help features
    - Call functions
    - Understand RGB colors

- **Part 2:** Playing with Blocks
    - Use IPython's help features
    - Assign variables
    - `for` loops (both iterator and `range` style)
    - `if` statements
    - indexing

- **Part 3:** Building Blocks
    - Reading files
    - Writing a function

## Overview

These lessons use [`ipythonblocks`][] to teach the basic usage of Python.
Making the most of `ipythonblocks` requires an understanding of RGB colors
so students should start with the "Color" notebook and spend a little time
practicing mixing their own colors. They may wish to leave this notebook
open as they move on so they can refer back to their color experiments.

The "Playing with Blocks" lesson goes into Python syntax, especially
indexing, `for` loops, and `if` statements. Finally, the "Building Blocks"
lesson goes into reading files and building functions.

[`ipythonblocks.py`][] is packaged here with the lesson so there's nothing
to install.

## Lesson Plans

*Note: These notebooks have some code pre-filled, such as imports.
Explain to the students that even though this code is already written in the
notebook, they must still explicitly execute those cells for them to have
any effect.*

### Color

Two functions are pre-imported: `show_color` and `embed_colorpicker`.

- Demonstrate using IPython's help features to see what these do, then explain
  [RGB colors][RGB] using `show_color` or other aids.

- Let them loose on the exercise, then show them how to use the [color picker
  tool][colorpicker] available via the `embed_colorpicker` function (requires an
  internet connection).

    *Hint: The RGB definitions of the colors for the exercise are visible in the
    raw Markdown cell.*

### Playing with Blocks

In this notebook the `BlockGrid` class has been imported for students.
The exercises are in the [play_with_blocks_exercises.md][exercises] file.

#### Variables

- Use IPython's help features to look at `BlockGrid`'s doc string.

- Demonstrate how to make a new `BlockGrid` object and assign it to a variable.
  This is a good chance to explain keyword arguments.

- Show how to display the grid using the interactive echo and the
  `BlockGrid.show()` method.

- Exercise 1

#### Basic Indexing

- Assign the first block of `grid` to a variable and change it's color,
  then display the grid.

        block = grid[0, 0]
        block.rgb = (0, 0, 0)
        grid

- Explain Python's zero based indexing, the coordinate system of the grid
  and that indices are `[row, column]`.

- Exercise 2

- Exercise 3
    - You can use `[-1, -1]` to get the lower-right block and explain
      Python's negative indexing.

- Exercise 4

#### Basic Loops

That's enough changing blocks one at a time, now for loops!

- Set the color of every block to something using a `for` loop:

        for block in grid:
            block.rgb = (0, 0, 0)

    This will probably be the first introduction of Pythonic indentation,
    so talk about that.

    Then demonstrate doing the same thing with the `.animate` attribute,
    which will show the changes as they happen:

        for block in grid.animate:
            block.rgb = (12, 123, 234)

- Exercise 5

#### Introducing If

Now to add logic so we can make some blocks different colors from others.

- Show an example `if` statement by looping over the grid but changing only
  one row. This will involve introducing the `block.row` attribute.

        for block in grid.animate:
            if block.row == 2:
                block.rgb = (0, 0, 0)

    A couple of new things are introduced here:

    - Using `==` for comparison vs. `=` for assignment.
      You might take this opportunity to introduce all of the comparison
      operators.
    - Indenting again for the `if` block.

    Also mention the `block.col` attribute.

- Exercise 6

- Now for `and`/`or`. Demo using `or` to change color of multiple columns
  with one loop through.

        for block in grid.animate:
            if block.col == 2 or block.col == 4:
                block.rgb = (50, 50, 50)

- Exercise 7

- Show the students that blocks have `.red`, `.green`, and `.blue` attributes
  they can use see the value of individual block color channels. (These can
  also be used to change the color values one at a time.)

- Exercise 8

[`ipythonblocks`]: https://github.com/jiffyclub/ipythonblocks
[`ipythonblocks.py`]: ./ipythonblocks.py
[RGB]: http://en.wikipedia.org/wiki/RGB_color_model
[colorpicker]: http://www.colorpicker.com
[exercises]: ./playing_with_blocks_exercises.md
