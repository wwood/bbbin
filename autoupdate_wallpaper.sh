#!/bin/bash

wget -o /dev/null -q -O $HOME/wallpaper.html http://www.lostsheep.com/lj.php 2>/dev/null
wget -o /dev/null -q -O $HOME/wallpaper `grep "img src" $HOME/wallpaper.html |head -n1 |sed -e "s/.*img src=.\(.*\). alt.*/\1/"` 2>/dev/null

# testing
#eog /tmp/wallpaper &

# reset background with a non-existent image
# if you don't use this, then next line will not update the wallpaper
gconftool -s -t string /desktop/gnome/background/picture_filename /usr/share/backgrounds/dummy.jpg
# load new desktop background image
gconftool -s -t string /desktop/gnome/background/picture_filename $HOME/wallpaper

