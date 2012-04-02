#!/bin/bash

# setup bashrc
if [ -f ~/.bashrc.bak ]; then
    cp ~/.bashrc.bak ~/.bashrc
fi
cp ~/.bashrc ~/.bashrc.bak
cp bashrc ~/.bashrc

# Install extra packages
apt-get update
apt-get install metacity -y
apt-get install gcc -y

# Fix to redraw borders
echo "@metacity --replace" >> /etc/xdg/lxsession/Lubuntu/autostart

# Install kernprof
wget http://pypi.python.org/packages/source/l/line_profiler/line_profiler-1.0b3.tar.gz#md5=63fc2a757192eb5e577559cfdff5b831
tar xvzf line_profiler-1.0b3.tar.gz
cd line_profiler-1.0b3
python setup.py install
cd ..
rm -rf line_profiler-1.0b3*

# Install faulthandler
wget http://pypi.python.org/packages/source/f/faulthandler/faulthandler-2.1.tar.gz#md5=5bbcb86e15a2b2583a97626f4b8c8760
tar xvzf faulthandler-2.1.tar.gz
cd faulthandler-2.1
python setup.py install
cd ..
rm -rf faulthandler-2.1*
