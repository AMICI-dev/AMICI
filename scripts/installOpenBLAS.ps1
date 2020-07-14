Write-Host 'script installOpenBLAS.ps1 started'
# $VerbosePreference = "Continue" # display verbose messages
New-Item -Path 'C:\BLAS' -ItemType Directory -Force # create directory
# Enforce stronger cryptography
[System.Net.ServicePointManager]::SecurityProtocol = [System.Net.SecurityProtocolType]::Tls12
# $uri = 'https://sourceforge.net/projects/openblas/files/v0.3.6/OpenBLAS-0.3.6-x64.zip/download'
$uri = 'https://github.com/xianyi/OpenBLAS/archive/v0.3.10.zip'
$output = 'C:\BLAS\v0.3.10.zip'
# Invoke-WebRequest $uri -OutFile $output
$webclient = New-Object System.Net.WebClient
$webclient.DownloadFile($uri,"$output")
Expand-Archive -Path 'C:\BLAS\v0.3.10.zip' -DestinationPath 'C:\BLAS\OpenBLAS-v0.3.10' -Force # expand zip file
#Set-Location "C:\BLAS\OpenBLAS-v0.3.10\OpenBLAS-0.3.10" # change directory
.\compileBLAS.cmd
#Copy-Item "C:\BLAS\OpenBLAS-v0.3.10\OpenBLAS-0.3.10\lib\openblas.lib" -Destination "C:\BLAS\lib"
#Copy-Item "C:\BLAS\OpenBLAS-v0.3.10\OpenBLAS-0.3.10\lib\openblas.dll" -Destination "C:\BLAS\bin"
# New-Item -Path Env:BLAS_LIBS -Value "/LIBPATH:C:\BLAS\\lib libopenblas.lib" -Force # create environment variable
# New-Item -Path Env:BLAS_CFLAGS -Value "/IC:\BLAS\OpenBLAS-v0.3.10\OpenBLAS-0.3.10\include" -Force # create environment variable
[System.Environment]::SetEnvironmentVariable("BLAS_LIBS", "/LIBPATH:C:\BLAS\lib openblas.lib", [System.EnvironmentVariableTarget]::Machine)
[System.Environment]::SetEnvironmentVariable("BLAS_LIBS", "/LIBPATH:C:\BLAS\lib openblas.lib", [System.EnvironmentVariableTarget]::User)
[System.Environment]::SetEnvironmentVariable("BLAS_LIBS", "/LIBPATH:C:\BLAS\lib openblas.lib", [System.EnvironmentVariableTarget]::Process)
[System.Environment]::SetEnvironmentVariable("BLAS_CFLAGS", "/IC:\BLAS\OpenBLAS-v0.3.10\OpenBLAS-0.3.10", [System.EnvironmentVariableTarget]::Machine)
[System.Environment]::SetEnvironmentVariable("BLAS_CFLAGS", "/IC:\BLAS\OpenBLAS-v0.3.10\OpenBLAS-0.3.10", [System.EnvironmentVariableTarget]::User)
[System.Environment]::SetEnvironmentVariable("BLAS_CFLAGS", "/IC:\BLAS\OpenBLAS-v0.3.10\OpenBLAS-0.3.10\include", [System.EnvironmentVariableTarget]::Process)
# Get-ChildItem 'C:\BLAS\OpenBLAS-v0.3.10' -Recurse # check for files
Get-Item -Path Env:BLAS_* # check environment variables
$VerbosePreference = "SilentlyContinue" # don't display verbose messages
Write-Host 'script installOpenBLAS.ps1 completed'