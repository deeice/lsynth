[Setup]
AppName=LSynth
AppVerName=LSynth 3.1
DefaultDirName={pf}\LSynth
DefaultGroupName=LSynth
Compression=lzma
LicenseFile="{app}\license.txt"

[Files]
Source: "{app}\bin\lsynthcp.exe"; DestDir: "{app}\bin"; Flags: ignoreversion 
Source: "{app}\bin\lsynth.mpd"; DestDir: "{app}\bin"; Flags: ignoreversion 
Source: "{app}\LSynthExampleParts.zip"; DestDir: "{app}"; Flags: ignoreversion 
Source: "{app}\Constraints.zip"; DestDir: "{app}"; Flags: ignoreversion 
Source: "{app}\readme.txt"; DestDir: "{app}"; Flags: ignoreversion 
Source: "{app}\license.txt"; DestDir: "{app}"; Flags: ignoreversion 

[Registry]
Root: HKCU; Subkey: "Software\LSynth"; ValueName: "InstallDir"; ValueType: String; ValueData: "{app}\bin"; 

[Run]
Filename: "{app}\readme.txt"; Description: "View the README file"; Flags: postinstall shellexec skipifsilent 
Filename: "{app}\bin\lsynthcp.exe"; Description: "{cm:LaunchProgram,LSynth}"; 

[Icons]
Name: "{group}\LSynth"; Filename: "{app}\bin\lsynthcp.exe"; 
Name: "{userdesktop}\LSynth"; Filename: "{app}\bin\lsynthcp.exe"; Tasks: "desktopicon"; 

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; 

[CustomMessages]
NameAndVersion=%1 version %2
AdditionalIcons=Additional icons:
CreateDesktopIcon=Create a &desktop icon
CreateQuickLaunchIcon=Create a &Quick Launch icon
ProgramOnTheWeb=%1 on the Web
UninstallProgram=Uninstall %1
LaunchProgram=Launch %1
AssocFileExtension=&Associate %1 with the %2 file extension
AssocingFileExtension=Associating %1 with the %2 file extension...
