Write-Host "Setting up Emscripten environment..." -ForegroundColor Green

if (-not (Test-Path "emsdk")) {
    Write-Host "Cloning emsdk..." -ForegroundColor Yellow
    git clone https://github.com/emscripten-core/emsdk.git
} else {
    Write-Host "emsdk already exists, skipping clone" -ForegroundColor Green
}

Set-Location "emsdk"
.\emsdk install latest
.\emsdk activate latest
.\emsdk_env.bat

Write-Host "Checking emcc version..." -ForegroundColor Green
emcc --version
if ($LASTEXITCODE -ne 0) {
    Write-Host "ERROR: emcc not found or not working properly" -ForegroundColor Red
    exit 1
}

Set-Location ".."
Write-Host "Current directory: $(Get-Location)"

if (-not (Test-Path "build")) {
    New-Item -ItemType Directory -Path "build" | Out-Null
}
Set-Location "build"

Write-Host "Running CMake configuration..." -ForegroundColor Green
emcmake cmake -G Ninja ..
if ($LASTEXITCODE -ne 0) {
    Write-Host "ERROR: CMake configuration failed" -ForegroundColor Red
    exit 1
}

Write-Host "Building project..." -ForegroundColor Green
emmake ninja
if ($LASTEXITCODE -ne 0) {
    Write-Host "ERROR: Build failed" -ForegroundColor Red
    exit 1
}

Write-Host "Checking for output files..." -ForegroundColor Green
$jsFiles = Get-ChildItem -Path "." -Filter "*.js" -ErrorAction SilentlyContinue
$wasmFiles = Get-ChildItem -Path "." -Filter "*.wasm" -ErrorAction SilentlyContinue

if ($jsFiles.Count -eq 0 -and $wasmFiles.Count -eq 0) {
    Write-Host "No output files found!" -ForegroundColor Red
} else {
    Write-Host "Build successful! Output files:" -ForegroundColor Green
    if ($jsFiles.Count -gt 0) {
        Write-Host "JavaScript files:" -ForegroundColor Yellow
        $jsFiles | ForEach-Object { Write-Host "  $($_.Name)" }
    }
    if ($wasmFiles.Count -gt 0) {
        Write-Host "WASM files:" -ForegroundColor Yellow
        $wasmFiles | ForEach-Object { Write-Host "  $($_.Name)" }
    }
}

Write-Host "Press any key to continue..." -ForegroundColor Gray
$null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")