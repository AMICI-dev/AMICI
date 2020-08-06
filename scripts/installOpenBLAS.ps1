Write-Host 'script installOpenBLAS.ps1 started'
New-Item -Path 'C:\BLAS' -ItemType Directory -Force # create directory
# Enforce stronger cryptography
[System.Net.ServicePointManager]::SecurityProtocol = [System.Net.SecurityProtocolType]::Tls12
$uri = 'https://github.com/xianyi/OpenBLAS/archive/v0.3.10.zip'
$output = 'C:\BLAS\v0.3.10.zip'
$webclient = New-Object System.Net.WebClient
$webclient.DownloadFile($uri,"$output")
Expand-Archive -Path 'C:\BLAS\v0.3.10.zip' -DestinationPath 'C:\BLAS\OpenBLAS-v0.3.10' -Force # expand zip file
C:\Users\travis\build\AMICI\scripts\compileBLAS.cmd
New-Item -Path 'C:\BLAS\lib' -ItemType Directory -Force # create directory
Copy-Item "C:\BLAS\OpenBLAS-v0.3.10\OpenBLAS-0.3.10\lib\Release\openblas.lib" -Destination "C:\BLAS\lib" -Recurse
New-Item -Path 'C:\BLAS\bin' -ItemType Directory -Force # create directory
Copy-Item "C:\BLAS\OpenBLAS-v0.3.10\OpenBLAS-0.3.10\lib\openblas.dll" -Destination "C:\BLAS\bin" -Recurse
Get-ChildItem -Path "C:\BLAS" -Include "openblas.lib" -Recurse # check for file
Get-ChildItem -Path "C:\BLAS" -Include "openblas.dll" -Recurse # check for file
[System.Environment]::SetEnvironmentVariable("BLAS_LIBS", "/LIBPATH:C:\BLAS\lib openblas.lib", [System.EnvironmentVariableTarget]::Machine)
[System.Environment]::SetEnvironmentVariable("BLAS_LIBS", "/LIBPATH:C:\BLAS\lib openblas.lib", [System.EnvironmentVariableTarget]::User)
[System.Environment]::SetEnvironmentVariable("BLAS_LIBS", "/LIBPATH:C:\BLAS\lib openblas.lib", [System.EnvironmentVariableTarget]::Process)
[System.Environment]::SetEnvironmentVariable("BLAS_CFLAGS", "-IC:\BLAS\OpenBLAS-v0.3.10\OpenBLAS-0.3.10", [System.EnvironmentVariableTarget]::Machine)
[System.Environment]::SetEnvironmentVariable("BLAS_CFLAGS", "-IC:\BLAS\OpenBLAS-v0.3.10\OpenBLAS-0.3.10", [System.EnvironmentVariableTarget]::User)
[System.Environment]::SetEnvironmentVariable("BLAS_CFLAGS", "-IC:\BLAS\OpenBLAS-v0.3.10\OpenBLAS-0.3.10", [System.EnvironmentVariableTarget]::Process)
Get-Item -Path Env:BLAS_* # check environment variables
$VerbosePreference = "SilentlyContinue" # don't display verbose messages
Write-Host 'script installOpenBLAS.ps1 completed'
