Write-Host 'script installOpenBLAS.ps1 started'
$version = '0.3.19'
New-Item -Path 'C:\BLAS' -ItemType Directory -Force # create directory
# Enforce stronger cryptography
[System.Net.ServicePointManager]::SecurityProtocol = [System.Net.SecurityProtocolType]::Tls12
$uri = "https://github.com/xianyi/OpenBLAS/archive/v$version.zip"
$output = "C:\BLAS\v$version.zip"
$webclient = New-Object System.Net.WebClient
$webclient.DownloadFile($uri,"$output")
Expand-Archive -Path "C:\BLAS\v$version.zip" -DestinationPath "C:\BLAS\OpenBLAS-$version" -Force # expand zip file
cmd /c "scripts\compileBLAS.cmd $version"
cmd /c dumpbin /DEPENDENTS "C:\BLAS\OpenBLAS-$version\OpenBLAS-$version\lib\openblas.dll"
New-Item -Path 'C:\BLAS\lib' -ItemType Directory -Force # create directory
Copy-Item "C:\BLAS\OpenBLAS-$version\OpenBLAS-$version\lib\Release\openblas.lib" -Destination "C:\BLAS\lib" -Recurse
New-Item -Path 'C:\BLAS\bin' -ItemType Directory -Force # create directory
Copy-Item "C:\BLAS\OpenBLAS-$version\OpenBLAS-$version\lib\openblas.dll" -Destination "C:\BLAS\bin" -Recurse
Get-ChildItem -Path "C:\BLAS" -Include "openblas.lib" -Recurse # check for file
Get-ChildItem -Path "C:\BLAS" -Include "openblas.dll" -Recurse # check for file
Get-Item -Path Env:BLAS_* # check environment variables
$VerbosePreference = "SilentlyContinue" # don't display verbose messages
Write-Host 'script installOpenBLAS.ps1 completed'
