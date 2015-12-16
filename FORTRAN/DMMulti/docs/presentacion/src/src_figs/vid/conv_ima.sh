#!/bin/bash

# Se convierte el video a imágener
#mplayer -vo jpeg -frames 50 -vf scale=360:416 $1 
#mplayer -vo jpeg -vf scale=216:250 $1 

# Se convierten los nombres del archivo
#rename 's/\d{5}(\d{3})\.jpg$/amo$1\.jpg/' *.jpg

# Ver bien cómo se manejan strings en bash
nom=$1
u=$((${#nom}-4))
ar=${nom:0:$u}
nombre="pic_${ar}%04d.jpg"
#ffmpeg -i "$1" -vf scale=w=iw/1.4:h=ih/1.4 $nombre
ffmpeg -i "$1" -vf scale=w=iw/1.2:h=ih/1.2 $nombre
