Get-ChildItem -Path C:\ -Include libopenblas.dll,python37.dll,MSVCP140.dll,KERNEL32.dll,VCRUNTIME140.dll,api-ms-win-crt-heap-l1-1-0.dll,api-ms-win-crt-stdio-l1-1-0.dll,api-ms-win-crt-string-l1-1-0.dll,api-ms-win-crt-runtime-l1-1-0.dll,api-ms-win-crt-math-l1-1-0.dll,api-ms-win-crt-time-l1-1-0.dll -File -Recurse -ErrorAction SilentlyContinue
Write-Host 'script installOpenBLAS.ps1 started'
$VerbosePreference = "Continue" # display verbose messages
New-Item -Path 'C:\BLAS' -ItemType Directory -Force # create directory
# Enforce stronger cryptography
[System.Net.ServicePointManager]::SecurityProtocol = [System.Net.SecurityProtocolType]::Tls12
$uri = 'https://sourceforge.net/projects/openblas/files/v0.3.6/OpenBLAS-0.3.6-x64.zip/download'
$output = 'C:\BLAS\OpenBLAS-0.3.6-x64.zip'
# Invoke-WebRequest $uri -OutFile $output
$webclient = New-Object System.Net.WebClient
$webclient.DownloadFile($uri,"$output")
Expand-Archive -Path 'C:\BLAS\OpenBLAS-0.3.6-x64.zip' -DestinationPath 'C:\BLAS\OpenBLAS-0.3.6-x64' -Force # expand zip file
New-Item -Path Env:BLAS_LIBS -Value "/LIBPATH:C:\BLAS\OpenBLAS-0.3.6-x64\lib libopenblas.lib" -Force # create environment variable
New-Item -Path Env:BLAS_CFLAGS -Value "/IC:\BLAS\OpenBLAS-0.3.6-x64\include" -Force # create environment variable
Get-ChildItem 'C:\BLAS\OpenBLAS-0.3.6-x64' -Recurse # check for files
Get-Item -Path Env:BLAS_* # check environment variables
$VerbosePreference = "SilentlyContinue" # don't display verbose messages
Write-Host 'script installOpenBLAS.ps1 completed'
