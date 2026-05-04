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
      message <- "Newton's Method cannot continue because f'(x_n) is too close to zero."
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
    tolerance = tolerance
  )
}

format_number <- function(value, digits = 8) {
  if (length(value) == 0 || is.na(value) || !is.finite(value)) {
    return("Not available")
  }

  formatC(value, format = "fg", digits = digits)
}

metric_box <- function(label, value) {
  div(
    class = "metric-box",
    span(class = "metric-label", label),
    span(class = "metric-value", value)
  )
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
        "family=Inter:wght@400;500;600;700&display=swap"
      )
    ),
    tags$style(HTML("
      :root {
        --ink: #18212f;
        --muted: #667085;
        --line: #d9e2ec;
        --panel: #ffffff;
        --paper: #f5f7f4;
        --accent: #157f73;
        --accent-dark: #0d5f57;
        --accent-soft: #e5f5f1;
        --warning: #a15c16;
        --danger: #a33131;
      }

      * {
        box-sizing: border-box;
      }

      body {
        margin: 0;
        color: var(--ink);
        background: var(--paper);
        font-family: Inter, Arial, sans-serif;
      }

      .container-fluid {
        max-width: 1240px;
        padding: 22px 24px 18px;
      }

      .app-shell {
        display: grid;
        gap: 18px;
      }

      .hero {
        background: #123c3a;
        color: #ffffff;
        border: 1px solid rgba(255,255,255,0.12);
        border-radius: 8px;
        padding: 22px 24px;
      }

      .hero h1 {
        margin: 0 0 7px;
        font-size: 30px;
        font-weight: 700;
      }

      .hero p {
        max-width: 760px;
        margin: 0;
        color: #d8ebe7;
        font-size: 15px;
        line-height: 1.55;
      }

      .layout {
        display: grid;
        grid-template-columns: 320px minmax(0, 1fr);
        gap: 18px;
        align-items: start;
      }

      .control-panel,
      .content-panel,
      .intro-panel {
        background: var(--panel);
        border: 1px solid var(--line);
        border-radius: 8px;
      }

      .control-panel {
        position: sticky;
        top: 18px;
        padding: 18px;
      }

      .control-panel h2,
      .content-panel h2,
      .intro-panel h2 {
        margin: 0 0 12px;
        font-size: 18px;
        font-weight: 700;
      }

      .form-group {
        margin-bottom: 13px;
      }

      label {
        font-size: 13px;
        font-weight: 650;
        color: #2e3a47;
      }

      .form-control {
        border-color: #cad6df;
        border-radius: 6px;
        box-shadow: none;
      }

      .form-control:focus {
        border-color: var(--accent);
        box-shadow: 0 0 0 3px rgba(21, 127, 115, 0.16);
      }

      .btn-primary {
        width: 100%;
        margin-top: 4px;
        background: var(--accent);
        border-color: var(--accent);
        border-radius: 6px;
        font-weight: 700;
      }

      .btn-primary:hover,
      .btn-primary:focus {
        background: var(--accent-dark);
        border-color: var(--accent-dark);
      }

      .input-note {
        margin: 12px 0 0;
        padding: 10px 11px;
        color: #405060;
        background: #f8faf7;
        border: 1px solid #e2e9e2;
        border-radius: 6px;
        font-size: 12px;
        line-height: 1.45;
      }

      .content-panel {
        padding: 0;
        overflow: hidden;
      }

      .nav-tabs {
        padding: 12px 14px 0;
        background: #fbfcfb;
        border-bottom: 1px solid var(--line);
      }

      .nav-tabs > li > a {
        color: #445265;
        border-radius: 6px 6px 0 0;
        font-weight: 650;
      }

      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:focus,
      .nav-tabs > li.active > a:hover {
        color: var(--accent-dark);
        border-color: var(--line);
        border-bottom-color: #ffffff;
      }

      .tab-content {
        padding: 18px;
        min-height: 560px;
      }

      .status-line {
        display: flex;
        gap: 10px;
        align-items: center;
        flex-wrap: wrap;
        margin-bottom: 14px;
      }

      .status-badge {
        display: inline-flex;
        align-items: center;
        min-height: 30px;
        padding: 5px 10px;
        border-radius: 999px;
        font-size: 13px;
        font-weight: 750;
      }

      .status-converged {
        color: #075c4f;
        background: var(--accent-soft);
      }

      .status-warning {
        color: var(--warning);
        background: #fff2df;
      }

      .status-error {
        color: var(--danger);
        background: #fdeaea;
      }

      .status-detail {
        color: var(--muted);
        font-size: 13px;
      }

      .metric-grid {
        display: grid;
        grid-template-columns: repeat(4, minmax(0, 1fr));
        gap: 10px;
        margin-bottom: 18px;
      }

      .metric-box {
        min-height: 88px;
        padding: 13px;
        background: #fbfcfb;
        border: 1px solid #e1e8df;
        border-radius: 8px;
      }

      .metric-label {
        display: block;
        margin-bottom: 8px;
        color: var(--muted);
        font-size: 12px;
        font-weight: 650;
        text-transform: uppercase;
      }

      .metric-value {
        display: block;
        color: var(--ink);
        font-size: 20px;
        font-weight: 750;
        word-break: break-word;
      }

      .method-note,
      .error-box {
        padding: 12px 14px;
        border-radius: 8px;
        line-height: 1.55;
      }

      .method-note {
        color: #435366;
        background: #f8faf7;
        border: 1px solid #e2e9e2;
      }

      .error-box {
        color: var(--danger);
        background: #fdeaea;
        border: 1px solid #facaca;
        font-weight: 650;
      }

      .intro-panel {
        padding: 18px;
      }

      .intro-grid {
        display: grid;
        grid-template-columns: 1fr 1fr;
        gap: 14px;
      }

      .intro-section {
        padding: 14px;
        background: #fbfcfb;
        border: 1px solid #e1e8df;
        border-radius: 8px;
      }

      .intro-section h3 {
        margin: 0 0 8px;
        font-size: 15px;
        font-weight: 750;
      }

      .intro-section p,
      .intro-section li {
        color: #465667;
        font-size: 14px;
        line-height: 1.55;
      }

      .intro-section ul {
        margin: 8px 0 0 18px;
        padding: 0;
      }

      .formula-strip {
        margin: 12px 0;
        padding: 16px;
        text-align: center;
        background: var(--accent-soft);
        border: 1px solid #c8e8df;
        border-radius: 8px;
      }

      .steps-wrap {
        display: grid;
        gap: 10px;
        max-height: 560px;
        overflow-y: auto;
        padding-right: 4px;
      }

      .step-box {
        padding: 13px 14px;
        background: #fbfcfb;
        border: 1px solid #e1e8df;
        border-left: 4px solid var(--accent);
        border-radius: 8px;
      }

      .step-box h4 {
        margin: 0 0 7px;
        font-size: 14px;
        font-weight: 750;
      }

      .step-box pre {
        margin: 0;
        padding: 0;
        color: #263342;
        background: transparent;
        border: 0;
        white-space: pre-wrap;
        font-size: 13px;
        line-height: 1.5;
      }

      .plot-caption {
        margin-top: 10px;
        color: var(--muted);
        font-size: 13px;
      }

      footer {
        color: var(--muted);
        text-align: center;
        font-size: 12px;
        padding: 8px 0 2px;
      }

      @media (max-width: 900px) {
        .container-fluid {
          padding: 14px;
        }

        .layout,
        .intro-grid {
          grid-template-columns: 1fr;
        }

        .control-panel {
          position: static;
        }

        .metric-grid {
          grid-template-columns: repeat(2, minmax(0, 1fr));
        }
      }

      @media (max-width: 560px) {
        .hero h1 {
          font-size: 24px;
        }

        .metric-grid {
          grid-template-columns: 1fr;
        }
      }
    "))
  ),
  div(
    class = "app-shell",
    div(
      class = "hero",
      h1("Newton's Method Solver"),
      p(
        "Explore how tangent-line approximations move from an initial ",
        "estimate toward a root of a nonlinear equation."
      )
    ),
    div(
      class = "layout",
      div(
        class = "control-panel",
        h2("Calculator"),
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
          "Initial estimate x0",
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
          "Use R syntax: 2*x, sin(x), cos(x), exp(x), log(x), sqrt(x), ",
          "and x^2. Leave the derivative blank to estimate it automatically."
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
              "Newton's Method uses the tangent line at x_n to choose the ",
              "next approximation. It is fast near a simple root, but it can ",
              "fail when the derivative is zero or the initial estimate is poor."
            )
          ),
          tabPanel(
            "Graph",
            plotOutput("function_plot", height = "500px"),
            div(
              class = "plot-caption",
              "The curve shows f(x), the horizontal line marks y = 0, ",
              "and the tangent line shows the Newton update from the latest iteration."
            )
          ),
          tabPanel(
            "Steps",
            uiOutput("steps")
          ),
          tabPanel(
            "Table",
            DT::dataTableOutput("table")
          ),
          tabPanel(
            "Introduction",
            div(
              class = "intro-panel",
              h2("Newton's Method"),
              p(
                "Newton's Method is an iterative technique for approximating ",
                "a root of f(x) = 0. Starting from x0, each step follows the ",
                "tangent line at the current point and uses where that tangent ",
                "crosses the x-axis as the next estimate."
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
                    "A practical stopping rule is based on the change between ",
                    "successive estimates: abs(x_{n+1} - x_n). Relative error ",
                    "uses that change divided by abs(x_{n+1})."
                  )
                ),
                div(
                  class = "intro-section",
                  h3("Input Rules"),
                  tags$ul(
                    tags$li("Write multiplication explicitly: 2*x, not 2x."),
                    tags$li("Use functions like sin(x), cos(x), exp(x), and log(x)."),
                    tags$li("Use powers with ^, such as x^6.")
                  )
                ),
                div(
                  class = "intro-section",
                  h3("Default Example"),
                  p(
                    "The app starts with f(x) = x^6 - x - 1 and x0 = 2, ",
                    "matching the reference example. It converges near 1.134724."
                  )
                )
              )
            )
          )
        )
      )
    ),
    tags$footer(
      "Numerical Analysis Project | Newton's Method"
    )
  )
)

server <- function(input, output, session) {
  calculation <- eventReactive(
    input$calculate,
    {
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
        metric_box("Root estimate", format_number(root_estimate, 10)),
        metric_box("f(root)", format_number(final_fx, 8)),
        metric_box("Approx. error", format_number(final_row$Approximate_Error, 8)),
        metric_box("Derivative mode", result$derivative_source)
      )
    )
  })

  output$table <- DT::renderDataTable({
    result <- calculation()
    validate(need(!isTRUE(result$error), result$message))

    display <- result$data
    names(display) <- c(
      "Iteration",
      "x_n",
      "f(x_n)",
      "f'(x_n)",
      "x_{n+1}",
      "Approximate Error",
      "Relative Error"
    )

    DT::datatable(
      display,
      rownames = FALSE,
      options = list(
        pageLength = 10,
        searching = FALSE,
        lengthChange = FALSE,
        scrollX = TRUE
      )
    ) |>
      DT::formatRound(
        columns = c("x_n", "f(x_n)", "f'(x_n)", "x_{n+1}", "Approximate Error", "Relative Error"),
        digits = 8
      )
  })

  output$steps <- renderUI({
    result <- calculation()

    if (isTRUE(result$error)) {
      return(div(class = "error-box", result$message))
    }

    rows <- result$data
    step_nodes <- lapply(seq_len(nrow(rows)), function(i) {
      row <- rows[i, ]

      if (!is.finite(row$x_next)) {
        text <- paste0(
          "x_", row$Iteration, " = ", format_number(row$x_n, 10), "\n",
          "f(x_", row$Iteration, ") = ", format_number(row$f_x_n, 10), "\n",
          "f'(x_", row$Iteration, ") = ", format_number(row$f_prime_x_n, 10), "\n",
          "The update cannot be computed because the derivative is too close to zero."
        )
      } else {
        text <- paste0(
          "x_", row$Iteration, " = ", format_number(row$x_n, 10), "\n",
          "f(x_", row$Iteration, ") = ", format_number(row$f_x_n, 10), "\n",
          "f'(x_", row$Iteration, ") = ", format_number(row$f_prime_x_n, 10), "\n",
          "x_", row$Iteration + 1, " = x_", row$Iteration,
          " - f(x_", row$Iteration, ") / f'(x_", row$Iteration, ")\n",
          "x_", row$Iteration + 1, " = ", format_number(row$x_n, 10),
          " - (", format_number(row$f_x_n, 10), ") / (",
          format_number(row$f_prime_x_n, 10), ")\n",
          "x_", row$Iteration + 1, " = ", format_number(row$x_next, 10), "\n",
          "Approximate error = ", format_number(row$Approximate_Error, 10)
        )
      }

      div(
        class = "step-box",
        h4(paste("Iteration", row$Iteration)),
        tags$pre(text)
      )
    })

    div(class = "steps-wrap", step_nodes)
  })

  output$function_plot <- renderPlot({
    result <- calculation()
    validate(need(!isTRUE(result$error), result$message))

    data <- result$data
    final_row <- data[nrow(data), ]
    root_estimate <- ifelse(
      is.finite(final_row$x_next),
      final_row$x_next,
      final_row$x_n
    )

    distance <- max(1, abs(input$x0 - root_estimate), abs(root_estimate) * 0.25)
    x_min <- min(input$x0, root_estimate) - distance
    x_max <- max(input$x0, root_estimate) + distance
    x_values <- seq(x_min, x_max, length.out = 500)

    y_values <- tryCatch(
      result$f(x_values),
      error = function(error) rep(NA_real_, length(x_values))
    )

    finite_points <- is.finite(y_values)
    validate(need(any(finite_points), "The function cannot be plotted on this range."))

    y_range <- range(y_values[finite_points], 0, na.rm = TRUE)
    y_padding <- diff(y_range) * 0.08

    if (!is.finite(y_padding) || y_padding == 0) {
      y_padding <- 1
    }

    plot(
      x_values,
      y_values,
      type = "l",
      lwd = 2,
      col = "#157f73",
      xlab = "x",
      ylab = "f(x)",
      main = "Newton's Method Graph",
      ylim = c(y_range[1] - y_padding, y_range[2] + y_padding)
    )
    grid(col = "#e5e9e2")
    abline(h = 0, col = "#606b78", lty = 2)

    points(
      input$x0,
      scalar_value(result$f, input$x0, "f(x0)"),
      pch = 19,
      col = "#a15c16",
      cex = 1.2
    )
    points(root_estimate, 0, pch = 19, col = "#123c3a", cex = 1.3)

    tangent_row <- data[nrow(data), ]

    if (is.finite(tangent_row$x_n) &&
        is.finite(tangent_row$f_x_n) &&
        is.finite(tangent_row$f_prime_x_n)) {
      tangent_y <- tangent_row$f_x_n +
        tangent_row$f_prime_x_n * (x_values - tangent_row$x_n)
      lines(x_values, tangent_y, col = "#c7524a", lwd = 2, lty = 3)
    }

    legend(
      "topright",
      legend = c("f(x)", "y = 0", "Initial estimate", "Root estimate", "Latest tangent"),
      col = c("#157f73", "#606b78", "#a15c16", "#123c3a", "#c7524a"),
      lty = c(1, 2, NA, NA, 3),
      pch = c(NA, NA, 19, 19, NA),
      bty = "n",
      cex = 0.9
    )
  })
}

shinyApp(ui, server)
