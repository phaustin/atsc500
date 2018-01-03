@echo off
rem this is a simple black and white prompt choice
PROMPT $P$Gcmd$_$G$S
set HOME=%USERPROFILE%
rem
rem add msys2 to your path
rem
set PATH=c:\msys64\usr\bin;%PATH%
rem
rem if mini36 isn't in your path ready, uncomment this lines
rem
rem set PATH=%userprofile%\mini36\scripts;%userprofile%\mini36;%userprofile%\mini36\Library\bin;%PATH%
set PWD=%~dp0
set h=%USERPROFILE%
@echo running %~dp0CmdInit.cmd
rem to make a shell that runs with this init file
rem to: cmd.exe /k "%UserProfile%\bin\CmdInit.cmd"

