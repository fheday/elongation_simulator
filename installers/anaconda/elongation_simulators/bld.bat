rem  load vstudio environment
rem set VSCMD_DEBUG=1
rem set VS_YEAR=2019
rem "C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\Common7\Tools\vsdevcmd"
rem if errorlevel 1 exit 1

rem build all simulators
rem msbuild "%SRC_DIR%"\elongation_simualtors.sln -p:Configuration=Release
rem if errorlevel 1 exit 1


cd installers\bin
rem if errorlevel 1 exit 1

copy vcruntime140.dll "%PREFIX%"\
if errorlevel 1 exit 1
copy ribosomesimulator.pyd "%PREFIX%"\
if errorlevel 1 exit 1
copy translation.pyd "%PREFIX%"\
if errorlevel 1 exit 1
mkdir "%PREFIX%"\concentrations
if errorlevel 1 exit 1
copy "%SRC_DIR%"\concentrations\__init__.py "%PREFIX%"\concentrations\
if errorlevel 1 exit 1
copy "%SRC_DIR%"\concentrations\Saccharomyces_cerevisiae.csv "%PREFIX%"\concentrations\
if errorlevel 1 exit 1

