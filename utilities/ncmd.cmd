@echo off
rem
rem  if ncmd.cmd is in a folder in %PATH% this will read in
rem  the cmd rcfile and cd to your home directory
rem
cmd.exe /k "%UserProfile%\bin\CmdInit.cmd && cd %USERPROFILE%"
