# load vstudio environment
#set VSCMD_DEBUG=1
#set VS_YEAR=2019
#"C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\Common7\Tools\vsdevcmd"
#if errorlevel 1 exit 1

#build all simulators
#msbuild "%SRC_DIR%"\elongation_simualtors.sln -p:Configuration=Release
#if errorlevel 1 exit 1


cd installers\bin
#if errorlevel 1 exit 1

copy ribosomesimulator.pyd "%PREFIX%"
if errorlevel 1 exit 1
copy translation.pyd "%PREFIX%"
if errorlevel 1 exit 1
#copy vcruntime140.dll "%PREFIX%"
#if errorlevel 1 exit 1

