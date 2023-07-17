@echo off
set oldPath=%cd%

cd /d %~dp0
cmake -S . -B build
cmake --open build

cd  /d %oldPath%