#!/bin/bash

#  Video encoding script
#  Author: Martin Hecht <mail@martin.hecht.de>
#  Date: 05/04/05

# default options:
FRAMERATE=10
BITRATE=8000
PLAYVIDEO=0
REMOVEPPMS=0
WIDTH=200
HEIGHT=200
#HEIGHT=288  

# catch filename if it is given as single parameter without '-o='...

AVIFILE=`echo $* | awk '{NTOT=NF; for (i = 1; i <= NTOT; i++) if (index($i, "-") != 1) print $i}'`

if [ "x$AVIFILE" != "x" ]
then
 echo "outputfile set to $AVIFILE"
fi


# function to print help message

# function to print help message

function printhelp
{
      cat <<- -EOHLP-

       mencode.sh  --  a script to build avi films 
                       it calls convert and mencoder

       mencode.sh [OPTIONS] 

       OPTIONS:
         
       -o=<filename>  : output file name for the video. 
                        default: output.avi

       -f=<framerate> : integer value of for the framerate, default: 25

       -b=<bitrate>   : integer value for the bitrate, default: 1800

       -w=<width>     : integer value for the image width, default: 384

       -h=<height>    : integer value for the image height, default: 288

       -p             : play the video after encoding, default: off

       -r             : remove jpg files after conversion, default: off

       -? or --help   : print this help message and exit

-EOHLP-
}


# find out all other command line parameters

PARAMETERS=`echo $* | awk '{NTOT=NF; for (i = 1; i <= NTOT; i++) if (index($i,"-") == 1) print $i}'`

  for p in $PARAMETERS
   do

    if echo $p | grep -e "-?" >/dev/null 2>&1
     then
       printhelp
       exit 0
     fi

    if echo $p | grep -e "--help" >/dev/null 2>&1
     then
       printhelp
       exit 0
     fi

    if echo $p | grep -e "-b" >/dev/null 2>&1
     then
       BITRATE=`echo $p | awk '{i=index($1, "="); print substr($1, i+1);}'`
       echo bitrate set to $BITRATE
     fi

    if echo $p | grep -e "-f" >/dev/null 2>&1
     then
       FRAMERATE=`echo $p | awk '{i=index($1, "="); print substr($1, i+1);}'`
       echo framerate set to $FRAMERATE
     fi

    if echo $p | grep -e "-o" >/dev/null 2>&1
     then
       AVIFILE=`echo $p | awk '{i=index($1, "="); print substr($1, i+1);}'`
       echo "outputfile set to $AVIFILE"        
     fi

    if echo $p | grep -e "-w" >/dev/null 2>&1
     then

       WIDTH=`echo $p | awk '{i=index($1, "="); print substr($1, i+1);}'`
       echo "width set to $WIDTH"        
     fi

    if echo $p | grep -e "-h" >/dev/null 2>&1
     then
       HEIGHT=`echo $p | awk '{i=index($1, "="); print substr($1, i+1);}'`
       echo "height set to $HEIGHT"        
     fi

    if echo $p | grep -e "-p" >/dev/null 2>&1
     then
      PLAYVIDEO=1
     fi

    if echo $p | grep -e "-r" >/dev/null 2>&1
     then
      REMOVEPPMS=1
     fi

  done


# now check if $AVIFILE is empty

if [ "x$AVIFILE" = "x" ]
 then
  AVIFILE="test.avi"
  echo "outputfile set to $AVIFILE"
 fi


# call mencoder with some encoding options

echo "calling now mencoder"

nice -19 mencoder mf://\*.jpg -mf w=${WIDTH}:h=${HEIGHT}:fps=${FRAMERATE}:type=jpg -ovc \
 lavc  -lavcopts vcodec=mpeg4:mbd=2:vbitrate=${BITRATE}:keyint=250 \
 -vf scale=${WIDTH}:${HEIGHT} -o $AVIFILE
# -vop scale=${WIDTH}:${HEIGHT} -o $AVIFILE



# and now play if the -p option was given

if [ $PLAYVIDEO -gt 0 ]
 then
  echo "encoding done, calling now mplayer"
  mplayer $AVIFILE
 fi

