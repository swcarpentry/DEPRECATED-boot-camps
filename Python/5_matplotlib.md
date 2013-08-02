
# Program Design as Exemplified in MatPlotLib

In the last session you saw how to use classes and objects to write software that can be used safely and easily by other developers. It is good practice to write software that can be used easily by other people. The best way to learn how to do this is to read and try to use other people's software. You will learn to improve your code by working out what they did right, and also by getting annoyed by the things that they did wrong.

MatPlotLib is a complete, object-orientated module for plotting graphs in Python. The code needed to draw graphs is obviously very complex. This complexity is hidden behind a set of classes / objects that are provided by MatPlotLib to provide a simple, matlab-like interface for users such as us.

First, load the MatPlotLib plotting library.

    $ ipython
    $ import matplotlib

The first thing to do when looking at some new code is to try to read the documentation...

    $ help(matplotlib)
    Help on package matplotlib:
    
    NAME
        matplotlib - This is an object-oriented plotting library.
    
    FILE
        /Users/chris/local/lib/python2.7/site-packages/matplotlib/__init__.py
    
    DESCRIPTION
        A procedural interface is provided by the companion pyplot module,
        which may be imported directly, e.g::
    
            from matplotlib.pyplot import *
    
        To include numpy functions too, use::
    
            from pylab import *
    
        or using ipython::
    
            ipython -pylab

This help is telling us that we want to use the "pyplot" companion module of matplotlib,

    $ from matplotlib import pyplot
    $ help(pyplot)
    Help on module matplotlib.pyplot in matplotlib:
    
    NAME
        matplotlib.pyplot - Provides a MATLAB-like plotting framework.
    
    FILE
        /Users/chris/local/lib/python2.7/site-packages/matplotlib/pyplot.py

You can see as I scroll down the help that it is very long... There are many classes and objects defined. As you can see, for big projects like this, simple code-level documentation is not enough. We need a simple quickstart guide. MatPlotLib provides this guide on their [website](http://matplotlib.org). This website provides a [tutorial](http://matplotlib.org/users/pyplot_tutorial.html), large numbers of [examples](http://matplotlib.org/examples/index.html) and a [gallery of graphs with code](http://matplotlib.org/gallery.html).

From the [tutorial](http://matplotlib.org/users/pyplot_tutorial.html), we can see that one of the main functions is "plot"

    $ help(pyplot.plot)
    plot(*args, **kwargs)
        Plot lines and/or markers to the
        :class:`~matplotlib.axes.Axes`.  *args* is a variable length
        argument, allowing for multiple *x*, *y* pairs with an
        optional format string.  For example, each of the following is
        legal::
    
            plot(x, y)         # plot x and y using default line style and color
            plot(x, y, 'bo')   # plot x and y using blue circle markers
            plot(y)            # plot y using x as index array 0..N-1
            plot(y, 'r+')      # ditto, but with red plusses


Let's now try to plot a simple y = x^2 graph using the "plot" command. We will use the example in the help and draw it with blue circles.

    $ x = []
    $ y = []
    $ for i in range(0,11):
    $     x.append(i)
    $     y.append(i*i)
    $
    $ pyplot.plot(x,y,'bo')

Strange - nothing is plotted. Looking back at the tutorial, we can see that we need to add use the "show" function.

    $ help(pyplot.show)
    Help on function show in module matplotlib.pyplot:
    
    show(*args, **kw)
        Display a figure.
    
        When running in ipython with its pylab mode, display all
        figures and return to the ipython prompt.

Now add "pyplot.show()"

    $ pyplot.show()

That should display the graph.

## Exercise

### Exercise 5a

Read through the [pyplot tutorial](http://matplotlib.org/users/pyplot_tutorial.html), typing in each of the examples. Use "help()" in ipython to understand each object and function that you are using.

Note, the tutorial imports the pyplot module using the command

    $ import matplotlib.pyplot as plt

This imports pyplot and renames it as "plt", e.g. so they can use the shorthand "plt.plot(x,y)" in place of "pyplot.plot(x,y)" etc.

## Object Orientated Frameworks

MatPlotLib is designed as a collection of interacting objects. Important objects are "Figure" and "Axes". The objects are designed and named to represent parts of a graph that a human can understand. The "Figure" class represents a graphical figure. The "Axes" object represents the axes of the graph. A good object-orientated design will involve objects that abstract human-understandable components of the program. For example, a molecular simulation program may contain Molecules and Atoms, while a crytography library may contain Keys and Signatures.

In a good object-orientated design, the functions and capabilities of the objects match their name, and would do what a user would expect would be reasonable. Lets use these objects to draw a graph...

    $ ipython
    $ from matplotlib import pyplot
    $ 
    $ figure = pyplot.figure()
    $ axes = figure.add_subplot(111)

This has returned the axes of the graph we will plot. We want to set the range and labels of the x- and y-axes. In most object orientated designs, functions that set variables are called "setSomething" or "set_something". Lets use ipython TAB to find all of the things we can set...

    $ axes.set[TAB]
    axes.set                       axes.set_clip_path             axes.set_transform
    axes.set_adjustable            axes.set_color_cycle           axes.set_url
    axes.set_agg_filter            axes.set_contains              axes.set_visible
    axes.set_alpha                 axes.set_cursor_props          axes.set_xbound
    axes.set_anchor                axes.set_figure                axes.set_xlabel
    axes.set_animated              axes.set_frame_on              axes.set_xlim
    axes.set_aspect                axes.set_gid                   axes.set_xmargin
    axes.set_autoscale_on          axes.set_label                 axes.set_xscale
    axes.set_autoscalex_on         axes.set_lod                   axes.set_xticklabels
    axes.set_autoscaley_on         axes.set_navigate              axes.set_xticks
    axes.set_axes                  axes.set_navigate_mode         axes.set_ybound
    axes.set_axes_locator          axes.set_picker                axes.set_ylabel
    axes.set_axis_bgcolor          axes.set_position              axes.set_ylim
    axes.set_axis_off              axes.set_rasterization_zorder  axes.set_ymargin
    axes.set_axis_on               axes.set_rasterized            axes.set_yscale
    axes.set_axisbelow             axes.set_snap                  axes.set_yticklabels
    axes.set_clip_box              axes.set_subplotspec           axes.set_yticks
    axes.set_clip_on               axes.set_title                 axes.set_zorder

Can you see what needs to be called? How about axes.set_xlim and axes.set_xlabel...

    $ help(axes.set_xlim)
    Help on method set_xlim in module matplotlib.axes:
    
    set_xlim(self, left=None, right=None, emit=True, auto=False, **kw) method of matplotlib.axes.AxesSubplot instance
        Call signature::
    
           set_xlim(self, *args, **kwargs):
    
        Set the data limits for the xaxis
    
    $ help(axes.set_xlabel)
    Help on method set_xlabel in module matplotlib.axes:
    
    set_xlabel(self, xlabel, fontdict=None, labelpad=None, **kwargs) method of matplotlib.axes.AxesSubplot instance
        Call signature::
      
          set_xlabel(xlabel, fontdict=None, labelpad=None, **kwargs)
     
        Set the label for the xaxis.
    
Ok - We will plot trigonometric functions, so lets set the axes...

    $ axes.set_xlim(0,360)
    $ axes.set_ylim(-1,1)
    $ axes.set_xlabel("Angle / degrees")
    $ axes.set_ylabel("Y")

Now lets generate some data for the graph. We will plot sine and cosine graphs. These functions are provided in the "math" module.

    $ import math
    $ x = []
    $ sin = []
    $ cos = []
    $ for i in range(0,360):
    $     x.append(i)
    $     angle_in_rads = 2 * math.pi * i / 360
    $     cos.append( math.cos(angle_in_rads) )
    $     sin.append( math.sin(angle_in_rads) )

We need to plot the data on the axes. Again, lets look for a plot function using ipython TAB...

    $ axes.p[TAB]
    axes.patch       axes.pcolor      axes.pick        axes.plot        axes.psd         
    axes.patches     axes.pcolorfast  axes.pickable    axes.plot_date   
    axes.pchanged    axes.pcolormesh  axes.pie         axes.properties  
    $ help(axes.plot)
    Help on method plot in module matplotlib.axes:
    
    plot(self, *args, **kwargs) method of matplotlib.axes.AxesSubplot instance
        Plot lines and/or markers to the
        :class:`~matplotlib.axes.Axes`.  *args* is a variable length
        argument, allowing for multiple *x*, *y* pairs with an
        optional format string.  For example, each of the following is
        legal::
    
            plot(x, y)         # plot x and y using default line style and color
            plot(x, y, 'bo')   # plot x and y using blue circle markers
            plot(y)            # plot y using x as index array 0..N-1
            plot(y, 'r+')      # ditto, but with red plusses

Lets plot the sine using a blue solid line, and cosine using a red solid line.

    $ axes.plot(x, sin, "b-")
    $ axes.plot(y, cos, "r-")

Now lets display the figure. Lets look for a function that can display the figure...

    $ figure.d[TAB]
    figure.delaxes          figure.dpi_scale_trans  figure.draw_artist      
    figure.dpi              figure.draw           

Nothing there... Perhaps MatPlotLib uses "show" instead of "display" or "draw"...

    $ figure.s[TAB]
    figure.savefig            figure.set_clip_on        figure.set_frameon        figure.set_transform
    figure.sca                figure.set_clip_path      figure.set_gid            figure.set_url
    figure.set                figure.set_contains       figure.set_label          figure.set_visible
    figure.set_agg_filter     figure.set_dpi            figure.set_lod            figure.set_zorder
    figure.set_alpha          figure.set_edgecolor      figure.set_picker         figure.show
    figure.set_animated       figure.set_facecolor      figure.set_rasterized     figure.subplotpars
    figure.set_axes           figure.set_figheight      figure.set_size_inches    figure.subplots_adjust
    figure.set_canvas         figure.set_figure         figure.set_snap           figure.suppressComposite
    figure.set_clip_box       figure.set_figwidth       figure.set_tight_layout   figure.suptitle

There it is...

    $ figure.show()

## Exercise

### Exercise 5b

Look through the functions of Axes. Try to draw the cosine and sine graphs using different axis styles and labels. For example look at axes.set_xticks and axes.grid. Try to guess what the function will do before you look at the help. Then, look at the help, e.g. help(axes.grid), and try to use the function. Does it do what you expected? Do you think that these are good names for the functions? What functions would you like to see on Axes, and what would you call them? Is the documentation clear? How would you improve it?

# [Previous](4_object_orientation.md) [Up](python_and_good_programming_practice.md)
