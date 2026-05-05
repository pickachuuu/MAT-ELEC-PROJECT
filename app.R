library(shiny)
library(DT)
library(numDeriv)

near_zero_derivative <- 1e-12

make_user_function <- function(expression_text) {
  expression_text <- trimws(expression_text)

  if (!nzchar(expression_text)) {
    stop("Expression is empty.", call. = FALSE)
  }

  parsed <- tryCatch(
    parse(text = expression_text),
    error = function(error) {
      stop(
        paste("R could not parse the expression:", conditionMessage(error)),
        call. = FALSE
      )
    }
  )

  function(x) {
    env <- new.env(parent = baseenv())
    env$x <- x
    value <- eval(parsed, envir = env)

    if (!is.numeric(value)) {
      stop("Expression must return a numeric value.", call. = FALSE)
    }

    if (length(value) == 1 && length(x) > 1) {
      value <- rep(value, length(x))
    }

    if (length(value) != length(x)) {
      stop(
        "Expression must return either one value or one value for each x.",
        call. = FALSE
      )
    }

    value
  }
}

scalar_value <- function(func, x, label) {
  value <- tryCatch(
    func(x),
    error = function(error) {
      stop(
        paste(label, "could not be evaluated:", conditionMessage(error)),
        call. = FALSE
      )
    }
  )

  if (length(value) != 1 || !is.numeric(value) || !is.finite(value)) {
    stop(
      paste(label, "must return one finite numeric value at the current x."),
      call. = FALSE
    )
  }

  as.numeric(value)
}

newton_method <- function(func_text, derivative_text, x0, tolerance, max_iter) {
  if (!is.numeric(x0) || length(x0) != 1 || !is.finite(x0)) {
    stop("Initial estimate x0 must be a finite number.", call. = FALSE)
  }

  if (!is.numeric(tolerance) || length(tolerance) != 1 ||
      !is.finite(tolerance) || tolerance <= 0) {
    stop("Tolerance must be a positive finite number.", call. = FALSE)
  }

  if (!is.numeric(max_iter) || length(max_iter) != 1 ||
      !is.finite(max_iter) || max_iter < 1) {
    stop("Maximum iterations must be at least 1.", call. = FALSE)
  }

  max_iter <- as.integer(max_iter)
  f <- make_user_function(func_text)
  derivative_source <- "Typed derivative"

  if (nzchar(trimws(derivative_text))) {
    df <- make_user_function(derivative_text)
  } else {
    derivative_source <- "Estimated derivative"
    df <- function(x) {
      numDeriv::grad(
        function(z) scalar_value(f, z, "f(x)"),
        x
      )
    }
  }

  rows <- vector("list", max_iter)
  status <- "Maximum iterations reached"
  message <- "The method stopped because it reached the iteration limit."
  x_current <- x0

  for (iteration in seq_len(max_iter)) {
    row_index <- iteration - 1
    fx <- scalar_value(f, x_current, "f(x)")
    dfx <- scalar_value(df, x_current, "f'(x)")

    if (abs(dfx) < near_zero_derivative) {
      status <- "Derivative near zero"
      message <- "Newton's Method cannot continue because the derivative at the current estimate is too close to zero."
      rows[[iteration]] <- data.frame(
        Iteration = row_index,
        x_n = x_current,
        f_x_n = fx,
        f_prime_x_n = dfx,
        x_next = NA_real_,
        Approximate_Error = NA_real_,
        Relative_Error = NA_real_,
        check.names = FALSE
      )
      break
    }

    x_next <- x_current - (fx / dfx)

    if (!is.finite(x_next)) {
      status <- "Non-finite result"
      message <- "Newton's Method produced a non-finite value."
      rows[[iteration]] <- data.frame(
        Iteration = row_index,
        x_n = x_current,
        f_x_n = fx,
        f_prime_x_n = dfx,
        x_next = x_next,
        Approximate_Error = NA_real_,
        Relative_Error = NA_real_,
        check.names = FALSE
      )
      break
    }

    approximate_error <- abs(x_next - x_current)
    relative_error <- approximate_error / max(abs(x_next), .Machine$double.eps)

    rows[[iteration]] <- data.frame(
      Iteration = row_index,
      x_n = x_current,
      f_x_n = fx,
      f_prime_x_n = dfx,
      x_next = x_next,
      Approximate_Error = approximate_error,
      Relative_Error = relative_error,
      check.names = FALSE
    )

    if (approximate_error <= tolerance) {
      status <- "Converged"
      message <- "The approximate error is within the requested tolerance."
      break
    }

    x_current <- x_next
  }

  data <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])

  list(
    data = data,
    status = status,
    message = message,
    f = f,
    derivative_source = derivative_source,
    tolerance = tolerance,
    x0 = x0
  )
}

format_number <- function(value, digits = 5) {
  if (length(value) == 0 || is.na(value) || !is.finite(value)) {
    return("Not available")
  }

  rounded <- round(value, digits)

  if (abs(rounded) >= 1e6 || (abs(value) < 1e-5 && value != 0)) {
    return(formatC(value, format = "e", digits = digits))
  }

  pretty <- formatC(rounded, format = "f", digits = digits)
  sub("\\.?0+$", "", pretty)
}

metric_box <- function(label, value) {
  div(
    class = "metric-box",
    span(class = "metric-label", label),
    span(class = "metric-value", value)
  )
}

presets <- list(
  custom    = list(label = "— Custom —"),
  default   = list(label = "x⁶ − x − 1  (default)",          func = "x^6 - x - 1",       deriv = "6*x^5 - 1",         x0 = 2),
  cosx      = list(label = "cos(x) − x",                      func = "cos(x) - x",        deriv = "-sin(x) - 1",        x0 = 1),
  cubic     = list(label = "x³ − 2x − 5  (Burden & Faires)", func = "x^3 - 2*x - 5",     deriv = "3*x^2 - 2",          x0 = 2),
  xexp      = list(label = "x − e^(−x)",                      func = "x - exp(-x)",       deriv = "1 + exp(-x)",        x0 = 0.5),
  sqrt2     = list(label = "x² − 2  (√2)",                    func = "x^2 - 2",           deriv = "2*x",                x0 = 1.5),
  fail_cyc  = list(label = "Failure: cycling",                func = "x^3 - 2*x + 2",     deriv = "3*x^2 - 2",          x0 = 0),
  fail_zero = list(label = "Failure: f'(x₀) = 0",            func = "x^3",                deriv = "3*x^2",              x0 = 0),
  fail_div  = list(label = "Failure: divergence (atan)",      func = "atan(x)",           deriv = "1/(1+x^2)",          x0 = 1.5)
)

preset_choices <- function() {
  setNames(names(presets), vapply(presets, function(p) p$label, character(1)))
}

ui <- fluidPage(
  tags$head(
    tags$link(
      rel = "preconnect",
      href = "https://fonts.googleapis.com"
    ),
    tags$link(
      rel = "preconnect",
      href = "https://fonts.gstatic.com",
      crossorigin = "anonymous"
    ),
    tags$link(
      rel = "stylesheet",
      href = paste0(
        "https://fonts.googleapis.com/css2?",
        "family=Caveat:wght@700&",
        "family=Patrick+Hand&",
        "family=Crimson+Pro:wght@600;700&",
        "family=Inter:wght@500;600&",
        "family=JetBrains+Mono:wght@500&display=swap"
      )
    ),
    tags$style(HTML("
      :root {
        --board-deep: #1f2a26;
        --board: #2c3a32;
        --board-light: #3a4d44;
        --chalk: #f5f1e8;
        --chalk-soft: #cfd0c4;
        --chalk-dim: rgba(245, 241, 232, 0.55);
        --chalk-faint: rgba(245, 241, 232, 0.18);
        --chalk-yellow: #f4d35e;
        --chalk-pink: #f4978e;
        --chalk-blue: #a8dadc;
        --chalk-green: #95d5b2;
        --paper-cream: #fdf6e3;
        --paper-line: rgba(168, 197, 227, 0.5);
        --paper-margin: #e63946;
        --ink-navy: #1d3557;
        --washi: rgba(244, 211, 94, 0.7);
        --board-shadow: 0 18px 40px rgba(0, 0, 0, 0.35);
        --paper-shadow: 0 6px 14px rgba(0, 0, 0, 0.28), 0 2px 4px rgba(0, 0, 0, 0.18);
      }

      * { box-sizing: border-box; }
      html, body { min-height: 100%; }

      body {
        margin: 0;
        color: var(--chalk);
        background-color: var(--board-deep);
        background-image:
          radial-gradient(ellipse at 18% 8%, rgba(245,241,232,0.06) 0%, transparent 55%),
          radial-gradient(ellipse at 82% 92%, rgba(245,241,232,0.05) 0%, transparent 60%),
          url(\"data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' width='220' height='220'><filter id='n'><feTurbulence type='fractalNoise' baseFrequency='0.85' numOctaves='2' stitchTiles='stitch'/><feColorMatrix values='0 0 0 0 0.96  0 0 0 0 0.94  0 0 0 0 0.91  0 0 0 0.05 0'/></filter><rect width='100%25' height='100%25' filter='url(%23n)'/></svg>\");
        background-size: auto, auto, 220px 220px;
        font-family: 'Patrick Hand', 'Segoe Print', cursive;
        font-size: 17px;
        line-height: 1.5;
      }

      body::before {
        content: '';
        position: fixed;
        inset: 0;
        pointer-events: none;
        background-image:
          repeating-linear-gradient(75deg, transparent 0 80px, rgba(245,241,232,0.018) 80px 81px),
          repeating-linear-gradient(105deg, transparent 0 110px, rgba(245,241,232,0.014) 110px 111px);
        z-index: 0;
      }

      .container-fluid {
        position: relative;
        max-width: 1280px;
        padding: 24px 28px;
        min-height: 100vh;
        z-index: 1;
      }

      .app-shell {
        display: grid;
        grid-template-rows: auto minmax(0, 1fr);
        gap: 18px;
        min-height: calc(100vh - 48px);
      }

      /* HERO */
      .hero {
        position: relative;
        overflow: hidden;
        padding: 26px 30px 22px;
        background: linear-gradient(180deg, rgba(245,241,232,0.025), transparent);
        border: 2px dashed var(--chalk-dim);
        border-radius: 10px;
      }

      .hero-content {
        display: grid;
        grid-template-columns: minmax(0, 1fr) auto;
        gap: 24px;
        align-items: start;
        position: relative;
        z-index: 1;
      }

      .hero-left { min-width: 0; }

      .hero-kicker {
        display: inline-flex;
        align-items: center;
        margin-bottom: 10px;
        padding: 4px 14px;
        color: var(--chalk-yellow);
        background: transparent;
        border: 1.5px dashed var(--chalk-yellow);
        border-radius: 999px;
        font-family: 'Patrick Hand', cursive;
        font-size: 14px;
        letter-spacing: 0.06em;
        text-transform: uppercase;
      }

      .hero h1, .hero-title {
        margin: 0 0 6px;
        font-family: 'Caveat', cursive;
        font-size: 60px;
        font-weight: 700;
        line-height: 1;
        color: var(--chalk);
        letter-spacing: 0.5px;
        position: relative;
        display: inline-block;
        text-shadow:
          1px 0 0 rgba(245,241,232,0.18),
          -1px 0 0 rgba(245,241,232,0.18),
          0 1px 0 rgba(245,241,232,0.12);
      }

      .hero-title::after {
        content: '';
        display: block;
        width: 78%;
        height: 9px;
        margin-top: 2px;
        background-image: url(\"data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 200 8' preserveAspectRatio='none'><path d='M2 5 Q 50 1, 100 4 T 198 5' fill='none' stroke='%23f5f1e8' stroke-width='2.4' stroke-linecap='round' opacity='0.85'/></svg>\");
        background-size: 100% 100%;
        background-repeat: no-repeat;
      }

      .hero p {
        max-width: 640px;
        margin: 8px 0 0;
        color: var(--chalk-soft);
        font-family: 'Patrick Hand', cursive;
        font-size: 18px;
        line-height: 1.45;
      }

      .hero-highlights {
        display: flex;
        flex-wrap: wrap;
        gap: 10px;
        justify-content: flex-end;
        max-width: 480px;
        margin-top: 4px;
      }

      .hero-chip {
        display: inline-flex;
        align-items: center;
        padding: 5px 13px;
        color: var(--chalk-blue);
        background: transparent;
        border: 1.5px dotted var(--chalk-blue);
        border-radius: 999px;
        font-family: 'Patrick Hand', cursive;
        font-size: 14px;
        letter-spacing: 0.02em;
      }

      .hero-chip.formula {
        color: var(--chalk-yellow);
        border-color: var(--chalk-yellow);
        font-family: 'JetBrains Mono', ui-monospace, monospace;
        font-size: 12px;
      }

      .hero-apple {
        position: absolute;
        right: 24px;
        bottom: 18px;
        width: 64px;
        height: 64px;
        opacity: 0.92;
        transform: rotate(-12deg);
        z-index: 0;
        pointer-events: none;
      }

      /* LAYOUT */
      .layout {
        display: grid;
        grid-template-columns: 320px minmax(0, 1fr);
        gap: 18px;
        align-items: stretch;
        min-height: 0;
      }

      /* SIDEBAR */
      .control-panel {
        position: sticky;
        top: 18px;
        align-self: start;
        padding: 22px 22px 20px;
        background: rgba(245,241,232,0.025);
        border: 2px dashed var(--chalk-dim);
        border-radius: 10px;
        box-shadow: var(--board-shadow);
      }

      .control-panel h2 {
        margin: 0 0 18px;
        font-family: 'Caveat', cursive;
        font-size: 34px;
        font-weight: 700;
        color: var(--chalk);
        position: relative;
        display: inline-block;
        line-height: 1;
      }

      .control-panel h2::after {
        content: '';
        display: block;
        width: 100%;
        height: 8px;
        margin-top: 2px;
        background-image: url(\"data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 200 8' preserveAspectRatio='none'><path d='M2 4 Q 50 1, 100 5 T 198 4' fill='none' stroke='%23f4d35e' stroke-width='2.6' stroke-linecap='round' opacity='0.95'/></svg>\");
        background-size: 100% 100%;
        background-repeat: no-repeat;
      }

      .form-group { margin-bottom: 14px; }

      .control-panel label, .control-panel .control-label {
        display: block;
        margin-bottom: 4px;
        font-family: 'Caveat', cursive;
        font-size: 22px;
        font-weight: 700;
        color: var(--chalk);
        letter-spacing: 0.3px;
      }

      .control-panel label sub, .control-panel label sup { font-family: 'Caveat', cursive; }

      .form-control {
        height: 42px;
        padding: 8px 12px;
        background: var(--board-light);
        color: var(--chalk);
        border: 1.5px solid var(--chalk-dim);
        border-radius: 6px;
        box-shadow: none;
        font-family: 'Inter', system-ui, sans-serif;
        font-weight: 500;
        font-size: 15px;
        transition: border-color 120ms, box-shadow 120ms;
      }

      .form-control:focus {
        outline: none;
        border-color: var(--chalk-yellow);
        box-shadow: 0 0 0 3px rgba(244, 211, 94, 0.22);
      }

      .form-control::placeholder { color: var(--chalk-soft); opacity: 0.6; }

      /* BUTTON */
      .btn-primary, .btn.btn-primary {
        display: inline-block;
        width: 100%;
        margin-top: 10px;
        padding: 8px 18px;
        min-height: 54px;
        background: var(--chalk-yellow);
        color: var(--board-deep) !important;
        border: 2.5px solid var(--board-deep);
        border-radius: 6px;
        font-family: 'Caveat', cursive;
        font-size: 28px;
        font-weight: 700;
        letter-spacing: 0.5px;
        line-height: 1;
        box-shadow: 4px 4px 0 0 rgba(245, 241, 232, 0.85);
        transition: transform 100ms, box-shadow 100ms;
        cursor: pointer;
        text-shadow: none;
      }

      .btn-primary:hover, .btn-primary:focus {
        background: #fce28a;
        color: var(--board-deep) !important;
        border-color: var(--board-deep);
        transform: translate(-1px, -1px);
        box-shadow: 5px 5px 0 0 rgba(245, 241, 232, 0.9);
        outline: none;
      }

      .btn-primary:active {
        transform: translate(2px, 2px);
        box-shadow: 1px 1px 0 0 rgba(245, 241, 232, 0.85);
      }

      /* INPUT NOTE — POST-IT */
      .input-note {
        position: relative;
        margin: 22px 4px 6px;
        padding: 18px 16px 14px;
        color: var(--ink-navy);
        background: var(--chalk-yellow);
        border-radius: 2px;
        font-family: 'Patrick Hand', cursive;
        font-size: 15px;
        line-height: 1.45;
        transform: rotate(-1.5deg);
        box-shadow: var(--paper-shadow);
      }

      .input-note::before {
        content: '';
        position: absolute;
        left: 50%;
        top: -10px;
        width: 72px;
        height: 20px;
        background: var(--washi);
        transform: translateX(-50%) rotate(-3deg);
        box-shadow: 0 1px 2px rgba(0,0,0,0.18);
      }

      /* CONTENT PANEL */
      .content-panel {
        padding: 0;
        overflow: hidden;
        display: flex;
        flex-direction: column;
        min-height: 0;
        background: rgba(245,241,232,0.025);
        border: 2px dashed var(--chalk-dim);
        border-radius: 10px;
        box-shadow: var(--board-shadow);
      }

      .content-panel > .tabbable {
        display: flex;
        flex: 1;
        flex-direction: column;
        min-height: 0;
      }

      /* TABS */
      .nav-tabs {
        display: flex;
        flex-wrap: wrap;
        gap: 18px;
        padding: 16px 24px 0;
        background: transparent;
        border: none;
        border-bottom: 1.5px dashed var(--chalk-faint);
        margin: 0;
      }

      .nav-tabs::before, .nav-tabs::after { display: none; }
      .nav-tabs > li { margin: 0; float: none; }

      .nav-tabs > li > a {
        position: relative;
        margin: 0;
        padding: 8px 4px 10px;
        color: var(--chalk-soft);
        background: transparent !important;
        border: none !important;
        border-radius: 0;
        font-family: 'Caveat', cursive;
        font-size: 24px;
        font-weight: 700;
        letter-spacing: 0.3px;
        line-height: 1;
      }

      .nav-tabs > li > a:hover {
        color: var(--chalk);
        background: transparent !important;
      }

      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:focus,
      .nav-tabs > li.active > a:hover {
        color: var(--chalk-yellow);
        background: transparent !important;
        border: none !important;
      }

      .nav-tabs > li.active > a::after {
        content: '';
        position: absolute;
        left: -2px;
        right: -2px;
        bottom: -3px;
        height: 7px;
        background-image: url(\"data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 200 8' preserveAspectRatio='none'><path d='M2 5 Q 50 1, 100 4 T 198 5' fill='none' stroke='%23f4d35e' stroke-width='2.5' stroke-linecap='round'/></svg>\");
        background-size: 100% 100%;
        background-repeat: no-repeat;
      }

      .tab-content {
        display: flex;
        flex: 1;
        min-height: 0;
        padding: 28px 26px 24px;
        background: transparent;
      }

      .tab-content > .tab-pane { width: 100%; }
      .tab-content > .active { display: flex; flex-direction: column; }

      /* STATUS LINE + STAMP */
      .status-line {
        display: flex;
        gap: 18px;
        align-items: center;
        flex-wrap: wrap;
        margin-bottom: 24px;
      }

      .status-badge {
        display: inline-flex;
        align-items: center;
        padding: 8px 18px;
        color: var(--paper-margin);
        background: transparent;
        border: 2.5px solid currentColor;
        border-radius: 4px;
        font-family: 'Patrick Hand', cursive;
        font-size: 18px;
        font-weight: 700;
        letter-spacing: 1.5px;
        text-transform: uppercase;
        transform: rotate(-3deg);
        box-shadow: inset 0 0 0 4px rgba(31, 42, 38, 0); /* spacer */
      }

      .status-badge::after {
        content: '';
        position: absolute;
        inset: 3px;
        border: 1.5px solid currentColor;
        border-radius: 3px;
        opacity: 0.55;
        pointer-events: none;
      }

      .status-badge { position: relative; }

      .status-converged { color: var(--paper-margin); }
      .status-warning { color: var(--chalk-yellow); }
      .status-error { color: var(--chalk-pink); }

      .status-detail {
        color: var(--chalk-soft);
        font-family: 'Patrick Hand', cursive;
        font-size: 17px;
      }

      /* METRIC CARDS — INDEX CARDS */
      .metric-grid {
        display: grid;
        grid-template-columns: repeat(4, minmax(0, 1fr));
        gap: 18px;
        margin-bottom: 24px;
      }

      .metric-box {
        position: relative;
        padding: 22px 18px 18px;
        min-height: 110px;
        background: var(--paper-cream);
        background-image: repeating-linear-gradient(
          180deg,
          transparent 0 26px,
          var(--paper-line) 26px 27px
        );
        background-position: 0 30px;
        border: none;
        border-radius: 2px;
        box-shadow: var(--paper-shadow);
        color: var(--ink-navy);
        overflow: visible;
      }

      .metric-box::before {
        content: '';
        position: absolute;
        width: 56px;
        height: 18px;
        top: -8px;
        background: var(--washi);
        box-shadow: 0 1px 2px rgba(0,0,0,0.18);
      }

      .metric-box:nth-child(1) { transform: rotate(-1deg); }
      .metric-box:nth-child(2) { transform: rotate(0.7deg); }
      .metric-box:nth-child(3) { transform: rotate(-0.6deg); }
      .metric-box:nth-child(4) { transform: rotate(1.1deg); }

      .metric-box:nth-child(odd)::before { left: 14px; transform: rotate(-6deg); }
      .metric-box:nth-child(even)::before { right: 14px; transform: rotate(5deg); }

      .metric-label {
        display: block;
        margin-bottom: 8px;
        color: var(--ink-navy);
        opacity: 0.7;
        font-family: 'Patrick Hand', cursive;
        font-size: 13px;
        font-weight: 400;
        text-transform: uppercase;
        letter-spacing: 0.12em;
      }

      .metric-value {
        display: block;
        color: var(--ink-navy);
        font-family: 'Crimson Pro', 'Times New Roman', serif;
        font-size: 26px;
        font-weight: 700;
        word-break: break-word;
        line-height: 1.15;
      }

      /* METHOD NOTE / ERROR — POST-IT */
      .method-note, .error-box {
        position: relative;
        padding: 22px 20px 18px;
        border-radius: 2px;
        font-family: 'Patrick Hand', cursive;
        font-size: 16px;
        line-height: 1.5;
        box-shadow: var(--paper-shadow);
      }

      .method-note {
        color: var(--ink-navy);
        background: var(--chalk-yellow);
        transform: rotate(-0.7deg);
        margin: 8px 6px 0;
      }

      .error-box {
        color: var(--ink-navy);
        background: var(--chalk-pink);
        transform: rotate(0.6deg);
        margin: 8px 6px 0;
      }

      .method-note::before, .error-box::before {
        content: '';
        position: absolute;
        left: 50%;
        top: -10px;
        width: 84px;
        height: 22px;
        background: var(--washi);
        transform: translateX(-50%) rotate(-2deg);
        box-shadow: 0 1px 2px rgba(0,0,0,0.18);
      }

      .method-note strong, .error-box strong {
        font-family: 'Caveat', cursive;
        font-size: 22px;
        font-weight: 700;
      }

      /* GRAPH PANE — PINNED PAPER */
      .paper-pin {
        position: relative;
        padding: 24px 22px 18px;
        background: var(--paper-cream);
        border-radius: 3px;
        box-shadow: var(--paper-shadow);
        transform: rotate(-0.4deg);
      }

      .paper-pin::before, .paper-pin::after {
        content: '';
        position: absolute;
        width: 86px;
        height: 22px;
        top: -10px;
        background: var(--washi);
        box-shadow: 0 1px 2px rgba(0,0,0,0.18);
      }

      .paper-pin::before { left: 30px; transform: rotate(-4deg); }
      .paper-pin::after { right: 30px; transform: rotate(3deg); }

      .paper-pin .shiny-plot-output, .paper-pin img { background: var(--paper-cream); }

      .plot-caption {
        margin-top: 14px;
        color: var(--chalk-soft);
        font-family: 'Patrick Hand', cursive;
        font-size: 15px;
        text-align: center;
      }

      /* STEPS — LINED PAPER CARDS */
      .steps-wrap {
        display: grid;
        gap: 16px;
        max-height: calc(100vh - 270px);
        overflow-y: auto;
        padding: 4px 8px 4px 4px;
      }

      .steps-wrap::-webkit-scrollbar { width: 8px; }
      .steps-wrap::-webkit-scrollbar-track { background: rgba(245,241,232,0.05); border-radius: 4px; }
      .steps-wrap::-webkit-scrollbar-thumb { background: var(--chalk-dim); border-radius: 4px; }

      .step-box {
        position: relative;
        padding: 16px 20px 14px 38px;
        background: var(--paper-cream);
        border-radius: 2px;
        box-shadow: var(--paper-shadow);
        background-image:
          linear-gradient(90deg, transparent 28px, rgba(230,57,70,0.55) 28px, rgba(230,57,70,0.55) 29px, transparent 29px),
          repeating-linear-gradient(180deg, transparent 0 26px, var(--paper-line) 26px 27px);
        background-position: 0 0, 0 32px;
        color: var(--ink-navy);
      }

      .step-box:nth-child(odd) { transform: rotate(-0.25deg); }
      .step-box:nth-child(even) { transform: rotate(0.25deg); }

      .step-box h4 {
        margin: 0 0 8px;
        font-family: 'Caveat', cursive;
        font-size: 26px;
        font-weight: 700;
        color: var(--ink-navy);
        line-height: 1;
      }

      .math-line {
        margin: 4px 0;
        color: var(--ink-navy);
        font-size: 15px;
        line-height: 1.5;
      }

      .math-line .MathJax, .math-line .MathJax_Display, .math-line mjx-container {
        color: var(--ink-navy) !important;
      }

      .math-note {
        margin-top: 8px;
        color: var(--ink-navy);
        opacity: 0.7;
        font-family: 'Patrick Hand', cursive;
        font-size: 14px;
        line-height: 1.45;
      }

      /* INTRO PANE */
      .intro-panel {
        padding: 0;
        background: transparent;
        border: none;
        box-shadow: none;
        color: var(--chalk);
      }

      .intro-panel h2 {
        margin: 0 0 12px;
        font-family: 'Caveat', cursive;
        font-size: 38px;
        color: var(--chalk);
        line-height: 1;
      }

      .intro-panel > p {
        font-family: 'Patrick Hand', cursive;
        font-size: 17px;
        color: var(--chalk-soft);
        line-height: 1.55;
        margin: 6px 0 0;
      }

      .formula-strip {
        margin: 20px 6px;
        padding: 24px 22px;
        text-align: center;
        background: var(--paper-cream);
        border-radius: 3px;
        box-shadow: var(--paper-shadow);
        transform: rotate(-0.3deg);
        position: relative;
        color: var(--ink-navy);
      }

      .formula-strip::before, .formula-strip::after {
        content: '';
        position: absolute;
        width: 64px;
        height: 18px;
        top: -8px;
        background: var(--washi);
        box-shadow: 0 1px 2px rgba(0,0,0,0.18);
      }
      .formula-strip::before { left: 40px; transform: rotate(-4deg); }
      .formula-strip::after { right: 40px; transform: rotate(3deg); }

      .formula-strip .MathJax,
      .formula-strip .MathJax_Display,
      .formula-strip mjx-container {
        color: var(--ink-navy) !important;
      }

      .intro-grid {
        display: grid;
        grid-template-columns: 1fr 1fr;
        gap: 20px;
        margin-top: 6px;
      }

      .intro-section {
        position: relative;
        padding: 22px 18px 16px;
        background: var(--paper-cream);
        border-radius: 2px;
        box-shadow: var(--paper-shadow);
        color: var(--ink-navy);
      }

      .intro-section:nth-child(1) { transform: rotate(-0.7deg); }
      .intro-section:nth-child(2) { transform: rotate(0.6deg); }
      .intro-section:nth-child(3) { transform: rotate(0.4deg); }
      .intro-section:nth-child(4) { transform: rotate(-0.5deg); }

      .intro-section::before {
        content: '';
        position: absolute;
        width: 54px;
        height: 16px;
        top: -7px;
        left: 18px;
        background: var(--washi);
        transform: rotate(-5deg);
        box-shadow: 0 1px 2px rgba(0,0,0,0.18);
      }

      .intro-section:nth-child(even)::before { left: auto; right: 18px; transform: rotate(5deg); }

      .intro-section h3 {
        margin: 0 0 8px;
        font-family: 'Caveat', cursive;
        font-size: 26px;
        font-weight: 700;
        color: var(--ink-navy);
        line-height: 1;
      }

      .intro-section p, .intro-section li {
        color: var(--ink-navy);
        font-family: 'Patrick Hand', cursive;
        font-size: 15px;
        line-height: 1.55;
      }

      .intro-section ul { margin: 8px 0 0 22px; padding: 0; }

      /* DT TABLE — CHALK TABLE */
      .dataTables_wrapper {
        font-family: 'Patrick Hand', cursive;
        color: var(--chalk);
        padding: 8px 4px 4px;
      }

      table.dataTable {
        background: transparent !important;
        color: var(--chalk) !important;
        border-collapse: collapse !important;
        width: 100% !important;
      }

      table.dataTable thead th {
        color: var(--chalk-yellow) !important;
        background: transparent !important;
        border-bottom: 2px solid var(--chalk-dim) !important;
        font-family: 'Caveat', cursive !important;
        font-size: 22px !important;
        font-weight: 700;
        letter-spacing: 0.3px;
        padding: 12px 14px !important;
        text-align: left;
      }

      table.dataTable tbody td {
        color: var(--chalk) !important;
        background: transparent !important;
        border: none !important;
        border-bottom: 1px dashed rgba(245, 241, 232, 0.22) !important;
        font-family: 'Patrick Hand', cursive !important;
        font-size: 17px !important;
        font-weight: 400 !important;
        padding: 11px 14px !important;
        letter-spacing: 0.2px;
      }

      table.dataTable tbody tr:hover td,
      table.dataTable.hover tbody tr:hover td {
        background: rgba(244, 211, 94, 0.10) !important;
        color: #ffffff !important;
      }

      table.dataTable.stripe tbody tr.odd td,
      table.dataTable.display tbody tr.odd > .sorting_1 {
        background: rgba(245,241,232,0.05) !important;
      }

      .dataTables_info, .dataTables_paginate, .dataTables_length, .dataTables_filter {
        color: var(--chalk-soft) !important;
        font-family: 'Patrick Hand', cursive !important;
        font-size: 15px !important;
        padding-top: 14px !important;
      }

      /* Pagination — three states: enabled, current, disabled.
         Selectors match DataTables' built-in specificity (.dataTables_wrapper ...). */
      .dataTables_wrapper .dataTables_paginate .paginate_button,
      .dataTables_wrapper .dataTables_paginate .paginate_button:link,
      .dataTables_wrapper .dataTables_paginate .paginate_button:visited {
        color: var(--chalk) !important;
        background: transparent !important;
        background-image: none !important;
        border: 1.5px solid var(--chalk-dim) !important;
        border-radius: 4px !important;
        margin: 0 3px;
        padding: 5px 14px !important;
        min-width: 36px;
        text-align: center;
        cursor: pointer;
        font-family: 'Patrick Hand', cursive !important;
        font-size: 16px !important;
        text-decoration: none !important;
        box-shadow: none !important;
        transition: all 120ms;
      }

      .dataTables_wrapper .dataTables_paginate .paginate_button:not(.disabled):not(.current):hover,
      .dataTables_wrapper .dataTables_paginate .paginate_button:not(.disabled):not(.current):active {
        color: var(--board-deep) !important;
        background: var(--chalk-yellow) !important;
        background-image: none !important;
        border-color: var(--chalk-yellow) !important;
        text-decoration: none !important;
      }

      .dataTables_wrapper .dataTables_paginate .paginate_button.current,
      .dataTables_wrapper .dataTables_paginate .paginate_button.current:link,
      .dataTables_wrapper .dataTables_paginate .paginate_button.current:visited,
      .dataTables_wrapper .dataTables_paginate .paginate_button.current:hover,
      .dataTables_wrapper .dataTables_paginate .paginate_button.current:active {
        color: var(--board-deep) !important;
        background: var(--chalk-yellow) !important;
        background-image: none !important;
        border: 2px solid var(--board-deep) !important;
        font-family: 'Caveat', cursive !important;
        font-size: 22px !important;
        font-weight: 700 !important;
        line-height: 1;
        box-shadow: 2px 2px 0 0 rgba(245, 241, 232, 0.55) !important;
        text-shadow: none !important;
      }

      .dataTables_wrapper .dataTables_paginate .paginate_button.disabled,
      .dataTables_wrapper .dataTables_paginate .paginate_button.disabled:hover,
      .dataTables_wrapper .dataTables_paginate .paginate_button.disabled:active {
        color: var(--chalk-soft) !important;
        background: transparent !important;
        background-image: none !important;
        border: 1.5px dashed rgba(245, 241, 232, 0.18) !important;
        opacity: 0.45 !important;
        cursor: not-allowed !important;
        box-shadow: none !important;
      }

      table.dataTable thead .sorting,
      table.dataTable thead .sorting_asc,
      table.dataTable thead .sorting_desc,
      table.dataTable thead .sorting_asc_disabled,
      table.dataTable thead .sorting_desc_disabled {
        position: relative;
        cursor: pointer;
        background-image: none !important;
        padding-right: 30px !important;
      }

      table.dataTable thead .sorting::before,
      table.dataTable thead .sorting::after,
      table.dataTable thead .sorting_asc::before,
      table.dataTable thead .sorting_asc::after,
      table.dataTable thead .sorting_desc::before,
      table.dataTable thead .sorting_desc::after,
      table.dataTable thead .sorting_asc_disabled::before,
      table.dataTable thead .sorting_asc_disabled::after,
      table.dataTable thead .sorting_desc_disabled::before,
      table.dataTable thead .sorting_desc_disabled::after {
        position: absolute !important;
        display: block !important;
        right: 10px;
        font-size: 11px !important;
        line-height: 1 !important;
        color: var(--chalk-soft) !important;
        opacity: 0.55 !important;
      }

      table.dataTable thead .sorting::before,
      table.dataTable thead .sorting_asc::before,
      table.dataTable thead .sorting_desc::before {
        content: \"\\25B2\";
        bottom: 54%;
      }

      table.dataTable thead .sorting::after,
      table.dataTable thead .sorting_asc::after,
      table.dataTable thead .sorting_desc::after {
        content: \"\\25BC\";
        top: 54%;
      }

      table.dataTable thead .sorting_asc::before {
        color: var(--chalk-yellow) !important;
        opacity: 1 !important;
      }

      table.dataTable thead .sorting_desc::after {
        color: var(--chalk-yellow) !important;
        opacity: 1 !important;
      }

      table.dataTable thead .sorting:hover::before,
      table.dataTable thead .sorting:hover::after {
        color: var(--chalk) !important;
        opacity: 0.85 !important;
      }

      /* hide DT's cloned thead inside scrollBody (kept for column-width sync) */
      .dataTables_scrollBody thead,
      .dataTables_scrollBody thead tr,
      .dataTables_scrollBody thead th {
        visibility: hidden !important;
      }
      .dataTables_scrollBody thead th::before,
      .dataTables_scrollBody thead th::after {
        display: none !important;
      }

      /* CHALKBOARD-WIDE TEXT (Introduction MathJax) */
      .intro-panel .MathJax, .intro-panel .MathJax_Display, .intro-panel mjx-container {
        color: var(--chalk) !important;
      }
      .formula-strip .MathJax, .formula-strip .MathJax_Display, .formula-strip mjx-container { color: var(--ink-navy) !important; }

      /* SHINY VALIDATION + ERROR — PINK POST-IT (matches .error-box) */
      .shiny-output-error-validation, .shiny-output-error {
        position: relative;
        margin: 8px 6px 0;
        padding: 22px 20px 18px;
        color: var(--ink-navy);
        background: var(--chalk-pink);
        border-radius: 2px;
        font-family: 'Patrick Hand', cursive;
        font-size: 16px;
        line-height: 1.5;
        transform: rotate(0.6deg);
        box-shadow: var(--paper-shadow);
      }

      .shiny-output-error-validation::before, .shiny-output-error::before {
        content: '';
        position: absolute;
        left: 50%;
        top: -10px;
        width: 84px;
        height: 22px;
        background: var(--washi);
        transform: translateX(-50%) rotate(-2deg);
        box-shadow: 0 1px 2px rgba(0,0,0,0.18);
      }

      /* MARGIN DOODLES */
      .doodle {
        position: fixed;
        pointer-events: none;
        opacity: 0.32;
        z-index: 0;
      }

      .doodle-parabola { top: 220px; right: 18px; width: 90px; transform: rotate(8deg); }
      .doodle-arrow { bottom: 90px; left: 24px; width: 76px; transform: rotate(-15deg); }

      /* COMPACT BUTTON (Play/Pause, Download CSV) */
      .btn-chalk-small, .btn.btn-chalk-small,
      a.btn-chalk-small.shiny-download-link {
        display: inline-flex !important;
        align-items: center;
        gap: 6px;
        padding: 6px 14px !important;
        min-height: 38px;
        background: var(--chalk-yellow) !important;
        color: var(--board-deep) !important;
        border: 2px solid var(--board-deep) !important;
        border-radius: 5px !important;
        font-family: 'Caveat', cursive !important;
        font-size: 20px !important;
        font-weight: 700 !important;
        letter-spacing: 0.5px;
        line-height: 1;
        box-shadow: 3px 3px 0 0 rgba(245, 241, 232, 0.85);
        cursor: pointer;
        text-decoration: none !important;
        transition: transform 100ms, box-shadow 100ms;
      }
      .btn-chalk-small:hover, .btn-chalk-small:focus,
      a.btn-chalk-small.shiny-download-link:hover {
        background: #fce28a !important;
        color: var(--board-deep) !important;
        transform: translate(-1px, -1px);
        box-shadow: 4px 4px 0 0 rgba(245, 241, 232, 0.9);
        outline: none;
        text-decoration: none !important;
      }
      .btn-chalk-small:active {
        transform: translate(1px, 1px);
        box-shadow: 1px 1px 0 0 rgba(245, 241, 232, 0.85);
      }
      .btn-chalk-small i, .btn-chalk-small .fa, .btn-chalk-small svg {
        color: var(--board-deep);
      }

      /* TABLE TOOLBAR (above DT) */
      .table-toolbar {
        display: flex;
        justify-content: flex-end;
        margin: 0 0 12px;
        padding: 0 4px;
      }

      /* ANIMATION CONTROLS (Graph tab) */
      .anim-controls {
        display: flex;
        gap: 16px;
        align-items: flex-end;
        margin-top: 14px;
        padding: 0 6px;
      }
      .anim-controls .form-group { flex: 1 1 auto; margin: 0; }
      .anim-controls .shiny-input-container { width: 100%; }
      .anim-controls > .btn-chalk-small { flex: 0 0 auto; margin-bottom: 4px; }

      /* CONVERGENCE CALLOUT (sticky note) */
      .convergence-callout {
        position: relative;
        margin: 18px 6px 0;
        padding: 22px 22px 18px;
        color: var(--ink-navy);
        background: var(--chalk-yellow);
        border-radius: 2px;
        font-family: 'Patrick Hand', cursive;
        font-size: 16px;
        line-height: 1.5;
        transform: rotate(-0.6deg);
        box-shadow: var(--paper-shadow);
      }
      .convergence-callout::before {
        content: '';
        position: absolute;
        left: 50%;
        top: -10px;
        width: 80px;
        height: 22px;
        background: var(--washi);
        transform: translateX(-50%) rotate(-2deg);
        box-shadow: 0 1px 2px rgba(0,0,0,0.18);
      }
      .convergence-callout strong {
        font-family: 'Caveat', cursive;
        font-size: 24px;
        font-weight: 700;
      }
      .convergence-order {
        display: inline-block;
        margin: 0 4px;
        padding: 2px 10px;
        background: rgba(29, 53, 87, 0.1);
        border: 1.5px solid var(--ink-navy);
        border-radius: 4px;
        font-family: 'Crimson Pro', serif;
        font-weight: 700;
        font-size: 18px;
      }

      /* IRS SLIDER OVERRIDES (sliderInput uses ionRangeSlider) */
      .anim-controls .irs--shiny .irs-bar { background: var(--chalk-yellow); border-top-color: var(--chalk-yellow); border-bottom-color: var(--chalk-yellow); }
      .anim-controls .irs--shiny .irs-line { background: rgba(245, 241, 232, 0.18); border-color: var(--chalk-dim); }
      .anim-controls .irs--shiny .irs-handle { background: var(--chalk-yellow); border: 2px solid var(--board-deep); box-shadow: 0 2px 4px rgba(0,0,0,0.3); }
      .anim-controls .irs--shiny .irs-handle > i { display: none; }
      .anim-controls .irs--shiny .irs-from,
      .anim-controls .irs--shiny .irs-to,
      .anim-controls .irs--shiny .irs-single { background: var(--chalk-yellow); color: var(--board-deep); font-family: 'Patrick Hand', cursive; font-size: 14px; }
      .anim-controls .irs--shiny .irs-from::before,
      .anim-controls .irs--shiny .irs-to::before,
      .anim-controls .irs--shiny .irs-single::before { border-top-color: var(--chalk-yellow); }
      .anim-controls .irs--shiny .irs-min,
      .anim-controls .irs--shiny .irs-max { color: var(--chalk-soft); background: transparent; font-family: 'Patrick Hand', cursive; font-size: 13px; }
      .anim-controls .irs--shiny .irs-grid-text { color: var(--chalk-soft); font-family: 'Patrick Hand', cursive; }
      .anim-controls .irs--shiny .irs-grid-pol { background: var(--chalk-soft); }
      .anim-controls label.control-label,
      .anim-controls .control-label { color: var(--chalk); font-family: 'Caveat', cursive; font-size: 22px; }

      /* RESPONSIVE */
      @media (max-width: 900px) {
        .container-fluid { padding: 16px; }
        .layout, .intro-grid { grid-template-columns: 1fr; }
        .hero-content { grid-template-columns: 1fr; }
        .hero-highlights { justify-content: flex-start; }
        .hero-apple { display: none; }
        .control-panel { position: static; }
        .app-shell, .layout { min-height: auto; }
        .metric-grid { grid-template-columns: repeat(2, minmax(0, 1fr)); }
        .doodle { display: none; }
      }

      @media (max-width: 560px) {
        .hero h1, .hero-title { font-size: 40px; }
        .metric-grid { grid-template-columns: 1fr; }
        .metric-box, .step-box, .input-note, .method-note,
        .error-box, .paper-pin, .formula-strip, .intro-section,
        .status-badge { transform: none; }
        .nav-tabs { gap: 12px; }
        .nav-tabs > li > a { font-size: 20px; }
      }
    ")),
    tags$script(HTML("
      $(document).on('keydown', '#func, #deriv, #x0, #tol, #max_iter', function(e) {
        if (e.key === 'Enter') { e.preventDefault(); $('#calculate').click(); }
      });
      $(document).on('shiny:connected', function() {
        Shiny.addCustomMessageHandler('clickCalculate', function(_) {
          $('#calculate').click();
        });
      });
    "))
  ),
  HTML("<svg xmlns='http://www.w3.org/2000/svg' style='position:absolute;width:0;height:0;overflow:hidden' aria-hidden='true'>
    <defs>
      <symbol id='sym-apple' viewBox='0 0 64 64'>
        <path d='M32 22 C 24 18, 14 24, 16 36 C 18 50, 28 56, 32 56 C 36 56, 46 50, 48 36 C 50 24, 40 18, 32 22 Z' fill='#f4978e' stroke='#f5f1e8' stroke-width='1.4' stroke-linejoin='round'/>
        <path d='M32 22 C 32 18, 30 12, 26 9' fill='none' stroke='#95d5b2' stroke-width='2.2' stroke-linecap='round'/>
        <path d='M32 22 C 36 18, 40 17, 44 16' fill='none' stroke='#95d5b2' stroke-width='2' stroke-linecap='round' opacity='0.85'/>
        <ellipse cx='27' cy='30' rx='3.5' ry='2.2' fill='#f5f1e8' opacity='0.55'/>
      </symbol>
      <symbol id='sym-parabola' viewBox='0 0 100 70'>
        <path d='M5 60 Q 50 -8, 95 60' fill='none' stroke='#f5f1e8' stroke-width='1.8' stroke-linecap='round' stroke-dasharray='3 4'/>
        <line x1='5' y1='50' x2='95' y2='50' stroke='#f5f1e8' stroke-width='1.2' opacity='0.55'/>
        <line x1='50' y1='5' x2='50' y2='65' stroke='#f5f1e8' stroke-width='1.2' opacity='0.55'/>
      </symbol>
      <symbol id='sym-arrow' viewBox='0 0 80 40'>
        <path d='M6 22 Q 30 4, 64 22' fill='none' stroke='#f5f1e8' stroke-width='1.8' stroke-linecap='round'/>
        <path d='M64 22 L 56 16 M64 22 L 58 28' fill='none' stroke='#f5f1e8' stroke-width='1.8' stroke-linecap='round'/>
      </symbol>
    </defs>
  </svg>"),
  HTML("<svg class='doodle doodle-parabola' viewBox='0 0 100 70' aria-hidden='true'><use href='#sym-parabola'/></svg>"),
  HTML("<svg class='doodle doodle-arrow' viewBox='0 0 80 40' aria-hidden='true'><use href='#sym-arrow'/></svg>"),
  div(
    class = "app-shell",
    div(
      class = "hero",
      div(
        class = "hero-content",
        div(
          class = "hero-left",
          span(class = "hero-kicker", "Numerical Analysis"),
          h1(class = "hero-title", "Newton's Method Solver"),
          p(
            "Explore how tangent-line approximations move from an initial ",
            "estimate toward a root of a nonlinear equation."
          )
        ),
        div(
          class = "hero-highlights",
          span(
            class = "hero-chip formula",
            HTML("x<sub>n+1</sub> = x<sub>n</sub> &minus; f(x<sub>n</sub>) / f'(x<sub>n</sub>)")
          ),
          span(class = "hero-chip", "Iteration table"),
          span(class = "hero-chip", "Tangent graph"),
          span(class = "hero-chip", "Auto derivative option")
        )
      ),
      HTML("<svg class='hero-apple' viewBox='0 0 64 64' aria-hidden='true'><use href='#sym-apple'/></svg>")
    ),
    div(
      class = "layout",
      div(
        class = "control-panel",
        h2("Calculator"),
        selectInput(
          "preset",
          "Examples",
          choices = preset_choices(),
          selected = "default"
        ),
        textInput(
          "func",
          "Function f(x)",
          value = "x^6 - x - 1"
        ),
        textInput(
          "deriv",
          "Derivative f'(x), optional",
          value = "6*x^5 - 1"
        ),
        numericInput(
          "x0",
          HTML("Initial estimate x<sub>0</sub>"),
          value = 2,
          step = 0.1
        ),
        numericInput(
          "tol",
          "Tolerance",
          value = 0.00001,
          min = 0.0000001,
          step = 0.00001
        ),
        numericInput(
          "max_iter",
          "Maximum iterations",
          value = 20,
          min = 1,
          step = 1
        ),
        actionButton(
          "calculate",
          "Calculate Root",
          class = "btn-primary"
        ),
        div(
          class = "input-note",
          HTML(
            "Use R syntax: 2*x, sin(x), cos(x), exp(x), log(x), sqrt(x), ",
            "and x<sup>2</sup>. Leave the derivative blank to estimate it automatically."
          )
        )
      ),
      div(
        class = "content-panel",
        tabsetPanel(
          tabPanel(
            "Results",
            uiOutput("summary"),
            div(
              class = "method-note",
              strong("Method note: "),
              HTML("Newton's Method uses the tangent line at x<sub>n</sub> to choose the "),
              "next approximation. It is fast near a simple root, but it can ",
              "fail when the derivative is zero or the initial estimate is poor."
            )
          ),
          tabPanel(
            "Graph",
            div(
              class = "paper-pin",
              plotOutput("function_plot", height = "500px", click = "plot_click")
            ),
            div(
              class = "anim-controls",
              sliderInput(
                "step_view",
                "Iteration",
                min = 0, max = 0, value = 0, step = 1,
                width = "100%"
              ),
              actionButton("play_toggle", "▶ Play", class = "btn-chalk-small")
            ),
            div(
              class = "plot-caption",
              HTML("Click anywhere on the plot to set <strong>x<sub>0</sub></strong>. "),
              "Drag the slider or press Play to step through iterations. ",
              "The dashed red line is the tangent at the current step."
            )
          ),
          tabPanel(
            "Convergence",
            div(class = "paper-pin", plotOutput("convergence_plot", height = "420px")),
            uiOutput("convergence_summary")
          ),
          tabPanel(
            "Steps",
            uiOutput("steps")
          ),
          tabPanel(
            "Table",
            div(
              class = "table-toolbar",
              downloadButton("download_csv", "Download CSV", class = "btn-chalk-small")
            ),
            DT::dataTableOutput("table")
          ),
          tabPanel(
            "Introduction",
            div(
              class = "intro-panel",
              h2("Newton's Method"),
              p(
                HTML(
                  "Newton's Method is an iterative technique for approximating ",
                  "a root of f(x) = 0. Starting from x<sub>0</sub>, each step follows the ",
                  "tangent line at the current point and uses where that tangent ",
                  "crosses the x-axis as the next estimate."
                )
              ),
              withMathJax(
                div(
                  class = "formula-strip",
                  "$$x_{n+1}=x_n-\\frac{f(x_n)}{f'(x_n)}, \\quad n \\ge 0$$"
                )
              ),
              div(
                class = "intro-grid",
                div(
                  class = "intro-section",
                  h3("Convergence"),
                  p(
                    "When f, f', and f'' are continuous near a root and ",
                    "f'(root) is not zero, Newton's Method can converge very ",
                    "quickly from a sufficiently close initial estimate."
                  )
                ),
                div(
                  class = "intro-section",
                  h3("Error Estimate"),
                  p(
                    HTML(
                      "A practical stopping rule is based on the change between ",
                      "successive estimates: |x<sub>n+1</sub> - x<sub>n</sub>|. Relative error ",
                      "uses that change divided by |x<sub>n+1</sub>|."
                    )
                  )
                ),
                div(
                  class = "intro-section",
                  h3("Input Rules"),
                  tags$ul(
                    tags$li("Write multiplication explicitly: 2*x, not 2x."),
                    tags$li("Use functions like sin(x), cos(x), exp(x), and log(x)."),
                    tags$li(HTML("Use powers with ^, such as x<sup>6</sup>."))
                  )
                ),
                div(
                  class = "intro-section",
                  h3("Default Example"),
                  p(
                    HTML(
                      "The app starts with f(x) = x<sup>6</sup> - x - 1 and x<sub>0</sub> = 2, ",
                      "matching the reference example. It converges near 1.13472."
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  calculation <- eventReactive(
    input$calculate,
    {
      validate(
        need(isTruthy(input$func), "Enter a function f(x)."),
        need(is.finite(input$x0), "Enter a finite number for the initial estimate x₀."),
        need(is.finite(input$tol) && input$tol > 0, "Tolerance must be a positive number."),
        need(is.finite(input$max_iter) && input$max_iter >= 1, "Maximum iterations must be at least 1.")
      )
      tryCatch(
        newton_method(
          func_text = input$func,
          derivative_text = input$deriv,
          x0 = input$x0,
          tolerance = input$tol,
          max_iter = input$max_iter
        ),
        error = function(error) {
          list(
            error = TRUE,
            status = "Invalid function/derivative",
            message = conditionMessage(error)
          )
        }
      )
    },
    ignoreInit = FALSE
  )

  # Helper: do current input values match the given preset?
  inputs_match_preset <- function(p) {
    if (is.null(p) || is.null(p$func)) return(FALSE)
    expected_deriv <- if (is.null(p$deriv)) "" else p$deriv
    isTruthy(input$func) && identical(input$func, p$func) &&
      identical(if (is.null(input$deriv)) "" else input$deriv, expected_deriv) &&
      isTruthy(input$x0) && is.finite(input$x0) &&
      abs(input$x0 - p$x0) < 1e-9
  }

  # When a preset is chosen, fill the form fields and re-run Calculate
  preset_applying <- reactiveVal(FALSE)

  observeEvent(input$preset, {
    if (input$preset == "custom") return()
    p <- presets[[input$preset]]
    if (is.null(p) || is.null(p$func)) return()
    if (inputs_match_preset(p)) {
      session$sendCustomMessage("clickCalculate", list())
      return()
    }
    preset_applying(TRUE)
    updateTextInput(session, "func", value = p$func)
    updateTextInput(session, "deriv", value = if (is.null(p$deriv)) "" else p$deriv)
    updateNumericInput(session, "x0", value = p$x0)
  }, ignoreInit = TRUE)

  # Reacts to any field edit. While a preset is being applied, don't flip to
  # Custom — instead wait for the inputs to converge to the preset's values,
  # then release the lock and trigger Calculate. After release, any edit flips
  # the dropdown to Custom.
  observeEvent(list(input$func, input$deriv, input$x0, input$tol, input$max_iter), {
    cur <- input$preset
    if (isTRUE(preset_applying())) {
      if (is.null(cur) || cur == "custom") {
        preset_applying(FALSE)
        return()
      }
      p <- presets[[cur]]
      if (inputs_match_preset(p)) {
        preset_applying(FALSE)
        session$sendCustomMessage("clickCalculate", list())
      }
      return()
    }
    if (!is.null(cur) && cur != "custom") {
      updateSelectInput(session, "preset", selected = "custom")
    }
  }, ignoreInit = TRUE)

  # Click on the function plot to set x0
  observeEvent(input$plot_click, {
    x_clicked <- round(input$plot_click$x, 4)
    if (!is.finite(x_clicked)) return()
    updateNumericInput(session, "x0", value = x_clicked)
    session$onFlushed(function() {
      session$sendCustomMessage("clickCalculate", list())
    }, once = TRUE)
  })

  output$summary <- renderUI({
    result <- calculation()

    if (isTRUE(result$error)) {
      return(div(class = "error-box", result$message))
    }

    data <- result$data
    final_row <- data[nrow(data), ]
    root_estimate <- ifelse(
      is.finite(final_row$x_next),
      final_row$x_next,
      final_row$x_n
    )
    final_fx <- tryCatch(
      scalar_value(result$f, root_estimate, "f(root)"),
      error = function(error) NA_real_
    )
    status_class <- switch(
      result$status,
      "Converged" = "status-converged",
      "Maximum iterations reached" = "status-warning",
      "Derivative near zero" = "status-error",
      "Non-finite result" = "status-error",
      "status-warning"
    )

    tagList(
      div(
        class = "status-line",
        span(class = paste("status-badge", status_class), result$status),
        span(class = "status-detail", result$message)
      ),
      div(
        class = "metric-grid",
        metric_box("Root estimate", format_number(root_estimate)),
        metric_box("f(root)", format_number(final_fx)),
        metric_box("Approx. error", format_number(final_row$Approximate_Error)),
        metric_box("Derivative mode", result$derivative_source)
      )
    )
  })

  output$table <- DT::renderDataTable({
    result <- calculation()
    validate(need(!isTRUE(result$error), result$message))

    display <- result$data
    display$Residual <- abs(display$f_x_n)
    display <- display[, c(
      "Iteration", "x_n", "f_x_n", "f_prime_x_n", "x_next",
      "Residual", "Approximate_Error", "Relative_Error"
    )]
    names(display) <- c(
      "Iteration",
      "x<sub>n</sub>",
      "f(x<sub>n</sub>)",
      "f'(x<sub>n</sub>)",
      "x<sub>n+1</sub>",
      "|f(x<sub>n</sub>)|",
      "Approximate Error",
      "Relative Error"
    )

    DT::datatable(
      display,
      class = "stripe hover compact cell-border",
      escape = FALSE,
      rownames = FALSE,
      options = list(
        pageLength = 10,
        searching = FALSE,
        lengthChange = FALSE,
        scrollX = TRUE
      )
    ) |>
      DT::formatRound(
        columns = 2:8,
        digits = 5
      )
  })

  output$download_csv <- downloadHandler(
    filename = function() {
      paste0("newton-iterations-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv")
    },
    content = function(file) {
      result <- calculation()
      req(!isTRUE(result$error))
      out <- result$data
      out$Residual <- abs(out$f_x_n)
      out <- out[, c(
        "Iteration", "x_n", "f_x_n", "f_prime_x_n", "x_next",
        "Residual", "Approximate_Error", "Relative_Error"
      )]
      write.csv(out, file, row.names = FALSE)
    }
  )

  output$steps <- renderUI({
    result <- calculation()

    if (isTRUE(result$error)) {
      return(div(class = "error-box", result$message))
    }

    rows <- result$data
    step_nodes <- lapply(seq_len(nrow(rows)), function(i) {
      row <- rows[i, ]
      current_index <- row$Iteration
      next_index <- current_index + 1

      if (!is.finite(row$x_next)) {
        step_body <- tagList(
          div(
            class = "math-line",
            paste0("\\(x_{", current_index, "} = ", format_number(row$x_n), "\\)")
          ),
          div(
            class = "math-line",
            paste0("\\(f(x_{", current_index, "}) = ", format_number(row$f_x_n), "\\)")
          ),
          div(
            class = "math-line",
            paste0("\\(\\left|f(x_{", current_index, "})\\right| = ", format_number(abs(row$f_x_n)), "\\)")
          ),
          div(
            class = "math-line",
            paste0("\\(f'(x_{", current_index, "}) = ", format_number(row$f_prime_x_n), "\\)")
          ),
          div(
            class = "math-note",
            "The update cannot be computed because the derivative is too close to zero."
          )
        )
      } else {
        step_body <- tagList(
          div(
            class = "math-line",
            paste0("\\(x_{", current_index, "} = ", format_number(row$x_n), "\\)")
          ),
          div(
            class = "math-line",
            paste0("\\(f(x_{", current_index, "}) = ", format_number(row$f_x_n), "\\)")
          ),
          div(
            class = "math-line",
            paste0("\\(\\left|f(x_{", current_index, "})\\right| = ", format_number(abs(row$f_x_n)), "\\)")
          ),
          div(
            class = "math-line",
            paste0("\\(f'(x_{", current_index, "}) = ", format_number(row$f_prime_x_n), "\\)")
          ),
          div(
            class = "math-line",
            paste0(
              "\\(x_{", next_index, "} = x_{", current_index,
              "} - \\frac{f(x_{", current_index, "})}{f'(x_{", current_index, "})}\\)"
            )
          ),
          div(
            class = "math-line",
            paste0(
              "\\(x_{", next_index, "} = ", format_number(row$x_n),
              " - \\frac{", format_number(row$f_x_n), "}{",
              format_number(row$f_prime_x_n), "}\\)"
            )
          ),
          div(
            class = "math-line",
            paste0("\\(x_{", next_index, "} = ", format_number(row$x_next), "\\)")
          ),
          div(
            class = "math-line",
            paste0("\\(\\left|x_{", next_index, "} - x_{", current_index, "}\\right| = ", format_number(row$Approximate_Error), "\\)")
          )
        )
      }

      div(
        class = "step-box",
        h4(paste("Iteration", current_index)),
        step_body
      )
    })

    withMathJax(div(class = "steps-wrap", step_nodes))
  })

  output$function_plot <- renderPlot({
    result <- calculation()
    validate(need(!isTRUE(result$error), result$message))

    data <- result$data
    n_total <- nrow(data) - 1L  # iterations 0..n_total
    step <- input$step_view
    if (is.null(step) || !is.finite(step)) step <- n_total
    step <- max(0L, min(as.integer(step), n_total))

    final_row <- data[nrow(data), ]
    root_estimate <- if (is.finite(final_row$x_next)) final_row$x_next else final_row$x_n

    x0 <- result$x0
    distance <- max(1, abs(x0 - root_estimate), abs(root_estimate) * 0.25)
    x_min <- min(x0, root_estimate) - distance
    x_max <- max(x0, root_estimate) + distance
    x_values <- seq(x_min, x_max, length.out = 500)

    y_values <- tryCatch(
      result$f(x_values),
      error = function(error) rep(NA_real_, length(x_values))
    )

    finite_points <- is.finite(y_values)
    validate(need(any(finite_points), "The function cannot be plotted on this range."))

    y_range <- range(y_values[finite_points], 0, na.rm = TRUE)
    y_padding <- diff(y_range) * 0.08
    if (!is.finite(y_padding) || y_padding == 0) y_padding <- 1

    op <- par(
      bg = "#fdf6e3", family = "serif",
      col.axis = "#1d3557", col.lab = "#1d3557", col.main = "#1d3557",
      fg = "#1d3557",
      cex.main = 1.4, cex.axis = 1, cex.lab = 1.1, font.main = 2
    )
    on.exit(par(op), add = TRUE)

    plot(
      x_values, y_values,
      type = "l", lwd = 2.6, col = "#1d3557",
      xlab = "x", ylab = "f(x)",
      main = paste0("Newton's Method  —  iteration ", step, " of ", n_total),
      ylim = c(y_range[1] - y_padding, y_range[2] + y_padding),
      panel.first = {
        usr <- par("usr")
        rect(usr[1], usr[3], usr[2], usr[4], col = "#fdf6e3", border = NA)
        grid(col = "#cfc6b8", lty = 1, lwd = 0.7)
      }
    )
    abline(h = 0, col = "#5b6b78", lty = 2, lwd = 1.4)

    # Iterates from row 1 (Iteration 0) up to and including current step
    visible_rows <- data[seq_len(step + 1L), , drop = FALSE]
    pt_x <- visible_rows$x_n
    pt_y <- vapply(pt_x, function(xv) tryCatch(scalar_value(result$f, xv, "f(x)"),
                                               error = function(e) NA_real_), numeric(1))
    finite_pts <- is.finite(pt_x) & is.finite(pt_y)
    if (any(finite_pts)) {
      points(pt_x[finite_pts], pt_y[finite_pts],
             pch = 19, col = "#1d3557", cex = 1.2)
      # Halo on the latest dot for emphasis
      latest <- which(finite_pts)[sum(finite_pts)]
      points(pt_x[latest], pt_y[latest],
             pch = 21, col = "#e63946", bg = "#1d3557", cex = 1.7, lwd = 2)
    }

    # Tangent at the current step
    cur_row <- data[step + 1L, ]
    if (is.finite(cur_row$x_n) && is.finite(cur_row$f_x_n) &&
        is.finite(cur_row$f_prime_x_n) && abs(cur_row$f_prime_x_n) > 1e-15) {
      tangent_y <- cur_row$f_x_n + cur_row$f_prime_x_n * (x_values - cur_row$x_n)
      lines(x_values, tangent_y, col = "#e63946", lwd = 2, lty = 3)
    }

    # Root marker only on the final frame
    if (step == n_total && is.finite(root_estimate)) {
      points(root_estimate, 0, pch = 19, col = "#e63946", cex = 1.5)
    }

    legend(
      "topright",
      legend = c("f(x)", "y = 0", "Iterate x_k", "Current iterate", "Tangent at x_k"),
      col = c("#1d3557", "#5b6b78", "#1d3557", "#e63946", "#e63946"),
      lty = c(1, 2, NA, NA, 3),
      pch = c(NA, NA, 19, 21, NA),
      pt.bg = c(NA, NA, NA, "#1d3557", NA),
      bg = "#fdf6e3", box.col = "#1d3557", text.col = "#1d3557",
      cex = 0.95
    )
  }, bg = "#fdf6e3")

  # Update the iteration slider whenever a new calculation finishes
  observeEvent(calculation(), {
    result <- calculation()
    if (isTRUE(result$error)) return()
    n_total <- nrow(result$data) - 1L
    updateSliderInput(session, "step_view", min = 0, max = n_total, value = n_total)
    playing(FALSE)
    updateActionButton(session, "play_toggle", label = "▶ Play")
  })

  # Animation engine
  playing <- reactiveVal(FALSE)

  observeEvent(input$play_toggle, {
    playing(!isTRUE(playing()))
    updateActionButton(
      session, "play_toggle",
      label = if (isTRUE(playing())) "⏸ Pause" else "▶ Play"
    )
  })

  observe({
    if (!isTRUE(playing())) return()
    invalidateLater(800, session)
    isolate({
      result <- calculation()
      if (isTRUE(result$error)) { playing(FALSE); return() }
      n_total <- nrow(result$data) - 1L
      cur <- input$step_view
      if (is.null(cur)) cur <- 0
      nxt <- if (cur >= n_total) 0L else as.integer(cur) + 1L
      updateSliderInput(session, "step_view", value = nxt)
    })
  })

  output$convergence_plot <- renderPlot({
    result <- calculation()
    validate(need(!isTRUE(result$error), result$message))

    data <- result$data
    final_row <- data[nrow(data), ]
    x_star <- if (is.finite(final_row$x_next)) final_row$x_next else final_row$x_n
    validate(need(is.finite(x_star), "No root estimate available."))

    errors <- abs(data$x_n - x_star)
    iter_idx <- data$Iteration
    keep <- is.finite(errors) & errors > 0
    errors <- errors[keep]
    iter_idx <- iter_idx[keep]
    validate(need(length(errors) >= 2, "Not enough iterations to plot convergence."))

    op <- par(
      bg = "#fdf6e3",
      family = "serif",
      col.axis = "#1d3557",
      col.lab = "#1d3557",
      col.main = "#1d3557",
      fg = "#1d3557",
      cex.main = 1.4, cex.axis = 1, cex.lab = 1.1,
      font.main = 2,
      mar = c(4.5, 5, 3.5, 1.5)
    )
    on.exit(par(op), add = TRUE)

    plot(
      iter_idx, errors,
      log = "y", type = "b", pch = 19, lwd = 2.4, cex = 1.2,
      col = "#1d3557",
      xlab = "Iteration n",
      ylab = expression(group("|", x[n] - x^"*", "|") ~ "(log scale)"),
      main = expression("Convergence:" ~ group("|", x[n] - x^"*", "|")),
      panel.first = {
        usr <- par("usr")
        rect(usr[1], 10^usr[3], usr[2], 10^usr[4], col = "#fdf6e3", border = NA)
        grid(col = "#cfc6b8", lty = 1, lwd = 0.7)
      }
    )

    show_guide <- length(errors) >= 3 && errors[2] < errors[1]
    if (show_guide) {
      # Quadratic guide: u_n = (2u_0 - u_1) + (u_1 - u_0) * 2^n,
      # i.e. a curve that's exact at n = 0, 1 and grows like e_n ~ e_{n-1}^2.
      u0 <- log(errors[1]); u1 <- log(errors[2])
      n_seq <- seq_along(errors) - 1
      log_ref <- (2 * u0 - u1) + (u1 - u0) * 2 ^ n_seq
      ref <- exp(pmax(log_ref, log(.Machine$double.eps)))
      lines(iter_idx, ref, col = "#e63946", lwd = 1.6, lty = 3)
    }
    legend(
      "topright",
      legend = if (show_guide)
        c(expression(group("|", x[n] - x^"*", "|")), "quadratic guide")
      else
        c(expression(group("|", x[n] - x^"*", "|"))),
      col = if (show_guide) c("#1d3557", "#e63946") else "#1d3557",
      lty = if (show_guide) c(1, 3) else 1,
      lwd = if (show_guide) c(2.4, 1.6) else 2.4,
      pch = if (show_guide) c(19, NA) else 19,
      bg = "#fdf6e3", box.col = "#1d3557", text.col = "#1d3557", cex = 0.95
    )
  }, bg = "#fdf6e3")

  output$convergence_summary <- renderUI({
    result <- calculation()
    if (isTRUE(result$error)) {
      return(div(class = "convergence-callout",
                 strong("Note: "), "the calculation did not produce a usable iteration sequence."))
    }
    data <- result$data
    final_row <- data[nrow(data), ]
    x_star <- if (is.finite(final_row$x_next)) final_row$x_next else final_row$x_n
    if (!is.finite(x_star)) {
      return(div(class = "convergence-callout",
                 strong("Note: "), "no finite root estimate is available for this run."))
    }
    errors <- abs(data$x_n - x_star)
    errors <- errors[is.finite(errors) & errors > 0]

    if (length(errors) < 3) {
      return(div(class = "convergence-callout",
                 strong("Not enough data: "),
                 "need at least three non-zero errors to estimate the order of convergence."))
    }

    n <- length(errors)
    ratios <- log(errors[(n-1):n]) / log(errors[(n-2):(n-1)])
    p_hat <- mean(ratios[is.finite(ratios)])
    p_hat <- max(0, min(p_hat, 4))
    p_text <- formatC(p_hat, format = "f", digits = 2)

    blurb <- if (p_hat >= 1.7) {
      "quadratic, as expected for Newton near a simple root."
    } else if (p_hat >= 1.2) {
      "between linear and quadratic — likely near a multiple root or far from regime."
    } else {
      "approximately linear — slow convergence (multiple root or weak conditions)."
    }

    div(class = "convergence-callout",
        strong("Estimated order: "),
        "p ≈ ", span(class = "convergence-order", p_text),
        " — ", blurb)
  })
}

shinyApp(ui, server)
