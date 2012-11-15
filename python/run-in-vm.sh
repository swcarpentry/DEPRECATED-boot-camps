#! /bin/bash
ipython notebook --pylab=inline --no-browser &
sleep 2
firefox http://127.0.0.1:8888 &
