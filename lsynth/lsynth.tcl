#!/bin/sh
#-*-tcl-*-
# the next line restarts using wish \
exec wish "$0" -- ${1+"$@"}

#
# Tk version of lsynth GUI
#

# wm withdraw .
wm title . "LSynth - LEGO Bendable Parts Synthesizer"

set ldrawdir $env(LDRAWDIR)
set modeldir $ldrawdir/models
set infile ""
set outfile ""

proc BrowseIF {} {
    global infile modeldir
    set infile [tk_getOpenFile -initialdir $modeldir]
}

proc BrowseOF {} {
    global outfile modeldir
    set outfile [tk_getOpenFile -initialdir $modeldir]
}

proc Synthesize {} {
    global infile outfile log
    if [catch {set n [exec bin/lsynthcp $infile $outfile]} input] {
      $log insert end $input\n
    } else {
      $log insert end n\n
    }
}

frame .f -borderwidth 10
pack .f -side top

frame .f.row0
pack .f.row0 -side top -anchor w
label .f.row0.li -text "Input File"
pack .f.row0.li -side left

frame .f.row1
pack .f.row1 -side top -fill y
entry .f.row1.inf -width 60 -relief sunken -textvariable infile
button .f.row1.browseInFile -text "Browse" -command {BrowseIF}
pack .f.row1.inf .f.row1.browseInFile -side left

frame .f.row2
pack .f.row2 -side top -anchor w
label .f.row2.lo -text "Output File"
pack .f.row2.lo -side left 

frame .f.row3
pack .f.row3 -side top -fill y
entry .f.row3.of -width 60 -relief sunken -textvariable outfile
button .f.row3.browseOutFile -text "Browse" -command {BrowseOF}
pack .f.row3.of .f.row3.browseOutFile -side left

frame .f.row4 -borderwidth 10
pack .f.row4 -side top -fill y
button .f.row4.synthesize -text "Synthesize" -command {Synthesize}
pack .f.row4.synthesize -side top

frame .f.row5
set log [text .f.row5.log -width 80 -height 20 \
    -borderwidth 2 -relief raised -setgrid true \
    -yscrollcommand {.f.row5.scroll set}]
scrollbar .f.row5.scroll -command {.f.row5.log yview}
pack .f.row5.scroll -side right -fill both -expand true
pack .f.row5 -side top -fill both -expand true

$log insert end "Welcome to LSynth, a bendable part synthesizer for LDRAW files."
$log insert end \n


