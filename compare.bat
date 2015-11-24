@echo off

rem === Set the python path to include sasview
for /D %%i in (..\sasview\build\lib*) do set SASVIEW=%%i
set PYTHONPATH=..\bumps;../periodictable;%SASVIEW%
rem echo PYTHONPATH=%PYTHONPATH%

rem === Set the pyopencl control parameters
rem set PYOPENCL_COMPILER_OUTPUT=1
rem set PYOPENCL_CTX=0

rem === Set OpenMP or not
rem set SAS_OPENMP=1

python -m sasmodels.compare %*
