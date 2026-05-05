#!/usr/bin/env bash
set -euo pipefail

APP="app.R"
PORT=3838
HOST_ADDRESS="127.0.0.1"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -a|--app) APP="$2"; shift 2 ;;
    -p|--port) PORT="$2"; shift 2 ;;
    -h|--host) HOST_ADDRESS="$2"; shift 2 ;;
    --help)
      echo "Usage: $0 [-a app.R] [-p 3838] [-h 127.0.0.1]"
      exit 0
      ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

if [[ -n "${R_SCRIPT:-}" && -x "${R_SCRIPT}" ]]; then
  RSCRIPT="$R_SCRIPT"
elif command -v Rscript >/dev/null 2>&1; then
  RSCRIPT="$(command -v Rscript)"
elif [[ -x "/Library/Frameworks/R.framework/Resources/bin/Rscript" ]]; then
  RSCRIPT="/Library/Frameworks/R.framework/Resources/bin/Rscript"
elif [[ -x "/opt/homebrew/bin/Rscript" ]]; then
  RSCRIPT="/opt/homebrew/bin/Rscript"
elif [[ -x "/usr/local/bin/Rscript" ]]; then
  RSCRIPT="/usr/local/bin/Rscript"
else
  echo "Rscript was not found. Install R (https://cran.r-project.org/bin/macosx/) or 'brew install r', or set R_SCRIPT to the Rscript path." >&2
  exit 1
fi

if [[ ! -f "$APP" ]]; then
  echo "App file not found: $APP" >&2
  exit 1
fi

RESOLVED_APP="$(cd "$(dirname "$APP")" && pwd)/$(basename "$APP")"

exec "$RSCRIPT" -e "renv::load('.'); if (!requireNamespace('shiny', quietly = TRUE)) { message('Restoring renv library (first run)...'); renv::restore(prompt = FALSE) }; shiny::runApp('${RESOLVED_APP}', host = '${HOST_ADDRESS}', port = ${PORT}, launch.browser = TRUE)"
