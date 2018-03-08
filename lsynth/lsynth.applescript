on run
  --prompt the user to choose an existing file
  set inputAlias to choose file
    try
      runLSynth(inputAlias)
    end try
end run


on open some_items
  repeat with inputAlias in some_items
    try
      runLSynth(inputAlias)
    end try
  end repeat
end open


to runLSynth(inputAlias)
set inputPath to POSIX path of inputAlias

--prompt for a file location
set defaultOutputName to name of (info for inputAlias)
set outputFile to choose file name default name ("Output " & defaultOutputName)
set outputPath to POSIX path of outputFile

--Find the command to execute. We have it stashed inside our script bundle 
set myPath to POSIX path of (path to me) as string

--pass it off to UNIX
-- (Unicode is preserved in the concatenation here)
-- the quoted form is documented at http://developer.apple.com/technotes/tn2002/tn2065.html#TNTAG4; 
-- but strangely not in the Standard Additions dictionary
set shellCommand to quoted form of myPath & "Contents/Resources/lsynthcp " & quoted form of inputPath & " " & quoted form of outputPath
do shell script shellCommand
set scriptOutput to the result

display dialog scriptOutput with title "LSynth Command Output" buttons {"OK"}
end runLSynth
