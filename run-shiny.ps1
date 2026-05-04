param(
  [string]$App = "references/sample.R",
  [int]$Port = 3838,
  [string]$HostAddress = "127.0.0.1"
)

$ErrorActionPreference = "Stop"

$ProjectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $ProjectRoot

$RscriptCandidates = @(
  $env:R_SCRIPT,
  "C:\Program Files\R\R-4.5.2\bin\Rscript.exe",
  "C:\Program Files\R\R-4.5.2\bin\x64\Rscript.exe"
) | Where-Object { $_ -and (Test-Path $_) }

if (-not $RscriptCandidates) {
  throw "Rscript.exe was not found. Install R for Windows or set the R_SCRIPT environment variable to Rscript.exe."
}

$ResolvedApp = (Resolve-Path $App).Path.Replace("\", "/")

& $RscriptCandidates[0] -e "renv::load('.'); shiny::runApp('$ResolvedApp', host = '$HostAddress', port = $Port, launch.browser = TRUE)"
