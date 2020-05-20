cd installers\bin
if errorlevel 1 exit 1

copy ribosomesimulator.pyd "%PREFIX%"
if errorlevel 1 exit 1
copy vcruntime140_1.dll "%PREFIX%"
if errorlevel 1 exit 1
copy vcruntime140.dll "%PREFIX%"
if errorlevel 1 exit 1

