cd installers\bin
rem if errorlevel 1 exit 1

copy vcruntime140_1.dll "%PREFIX%"
if errorlevel 1 exit 1
copy ribosomesimulator.pyd "%PREFIX%"
if errorlevel 1 exit 1
copy translation.pyd "%PREFIX%"\
if errorlevel 1 exit 1
mkdir "%PREFIX%"\concentrations
if errorlevel 1 exit 1
copy "%SRC_DIR%"\concentrations\__init__.py "%PREFIX%"\concentrations\
if errorlevel 1 exit 1
copy "%SRC_DIR%"\concentrations\Saccharomyces_cerevisiae.csv "%PREFIX%"\concentrations\
if errorlevel 1 exit 1

