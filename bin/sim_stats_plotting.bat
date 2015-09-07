@echo off

python %~dp0\sim_stats_plotting.py %*

IF ERRORLEVEL 1 pause
