$ErrorActionPreference = "Stop"

$projectRoot = Split-Path -Parent $PSScriptRoot
$python = "C:\Users\dell\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe"
$prepare = Join-Path $PSScriptRoot "prepare_zenodo_v1_2.py"
$publisher = Join-Path $PSScriptRoot "publish_zenodo_v1_2.py"
$log = Join-Path $projectRoot "public_repository\zenodo_v1.2.0_publish.log"

Set-Location -LiteralPath $projectRoot
Write-Host "Preparing Zenodo v1.2.0 archive..." -ForegroundColor Cyan
& $python $prepare
if ($LASTEXITCODE -ne 0) {
    throw "Archive preparation failed."
}

Write-Host ""
Write-Host "Paste the Zenodo token and press Enter." -ForegroundColor Yellow
Write-Host "No characters or asterisks will appear while you type or paste. This is normal."
$secureToken = Read-Host "Zenodo token" -AsSecureString
$tokenPointer = [Runtime.InteropServices.Marshal]::SecureStringToBSTR($secureToken)

try {
    $env:ZENODO_TOKEN = [Runtime.InteropServices.Marshal]::PtrToStringBSTR($tokenPointer)
    if ([string]::IsNullOrWhiteSpace($env:ZENODO_TOKEN)) {
        throw "No token was entered."
    }
    Write-Host "Token received. Creating and uploading Zenodo v1.2.0..." -ForegroundColor Green
    $previousPreference = $ErrorActionPreference
    $ErrorActionPreference = "Continue"
    & $python $publisher 2>&1 |
        ForEach-Object { $_.ToString() } |
        Tee-Object -FilePath $log
    $ErrorActionPreference = $previousPreference
    if ($LASTEXITCODE -ne 0) {
        throw "Zenodo publisher exited with code $LASTEXITCODE. See $log"
    }
    Write-Host "Zenodo v1.2.0 publication completed." -ForegroundColor Green
}
catch {
    $_ | Out-String | Tee-Object -FilePath $log -Append
    Write-Host "Publication failed. Keep this window open and ask Codex to inspect the log." -ForegroundColor Red
}
finally {
    $ErrorActionPreference = "Continue"
    [Runtime.InteropServices.Marshal]::ZeroFreeBSTR($tokenPointer)
    Remove-Item Env:ZENODO_TOKEN -ErrorAction SilentlyContinue
}
