#!/bin/sh

osacompile -o LSynth.app lsynth.applescript

cp lsynth/lsynthcp LSynth.app/Contents/Resources
cp lsynth/lsynth.mpd LSynth.app/Contents/Resources


