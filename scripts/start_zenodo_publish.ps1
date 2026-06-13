$ErrorActionPreference = "Stop"

$projectRoot = Split-Path -Parent $PSScriptRoot
$python = "C:\Users\dell\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe"
$publisher = Join-Path $PSScriptRoot "publish_zenodo.py"
$log = Join-Path $projectRoot "public_repository\zenodo_publish.log"

Set-Location -LiteralPath $projectRoot

Write-Host "Zenodo publication for TIMP1 research compendium" -ForegroundColor Cyan
Write-Host "The token will not be displayed or written to disk."
$secureToken = Read-Host "Paste Zenodo token" -AsSecureString
$tokenPointer = [Runtime.InteropServices.Marshal]::SecureStringToBSTR($secureToken)

try {
    $env:ZENODO_TOKEN = [Runtime.InteropServices.Marshal]::PtrToStringBSTR($tokenPointer)
    Write-Host "Starting Zenodo upload. The 230 MB archive may take several minutes..." -ForegroundColor Yellow
    # Windows PowerShell wraps native stderr as ErrorRecord objects. Do not let
    # those records terminate the launcher before Python prints its traceback.
    $previousErrorActionPreference = $ErrorActionPreference
    $ErrorActionPreference = "Continue"
    & $python $publisher 2>&1 |
        ForEach-Object { $_.ToString() } |
        Tee-Object -FilePath $log
    $ErrorActionPreference = $previousErrorActionPreference
    if ($LASTEXITCODE -ne 0) {
        throw "Zenodo publisher exited with code $LASTEXITCODE. See $log"
    }
    Write-Host ""
    Write-Host "Zenodo publication completed." -ForegroundColor Green
}
catch {
    $_ | Out-String | Tee-Object -FilePath $log -Append
    Write-Host ""
    Write-Host "Publication failed. Keep this window open and tell Codex to inspect the log." -ForegroundColor Red
}
finally {
    $ErrorActionPreference = "Continue"
    [Runtime.InteropServices.Marshal]::ZeroFreeBSTR($tokenPointer)
    Remove-Item Env:ZENODO_TOKEN -ErrorAction SilentlyContinue
}
