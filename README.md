# MAT-ELEC-PROJECT

Newton's Method Shiny app for a numerical analysis project.

## Run

Windows (PowerShell):

```powershell
.\run-shiny.ps1 -App "app.R"
```

macOS / Linux:

```bash
./run-shiny.sh --app app.R
```

The project uses `renv` for R package management. If packages need to be restored on another machine, run:

```r
renv::restore()
```
