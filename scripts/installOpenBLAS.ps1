Write-Host 'script installOpenBLAS.ps1 started'
$VerbosePreference = "Continue" # display verbose messages
New-Item -Path 'C:\BLAS' -ItemType Directory -Force # create directory
# Set-Location 'C:\BLAS'
# Copy-Item -Path "$BASE_DIR\ThirdParty\BLAS\OpenBLAS-0.3.6-x64.zip" -Destination "C:\BLAS" -Force
# Enforce stronger cryptography
# [Net.ServicePointManager]::SecurityProtocol = [Net.ServicePointManager]::SecurityProtocol -bor [Net.SecurityProtocolType]::Tls12
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
