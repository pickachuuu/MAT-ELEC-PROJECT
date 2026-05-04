library(shiny)
library(DT)
library(numDeriv)

central_difference <- function(f, x, h) {
  (f(x + h) - f(x - h)) / (2 * h)
}

# Define UI
ui <- fluidPage(
  tags$head(
    HTML(
      "<link href='https://fonts.googleapis.com/css2?family=Roboto:ital,wght@0,100;0,300;0,400;0,500;0,700;0,900;1,100;1,300;1,400;1,500;1,700;1,900&display=swap' rel='stylesheet'>"
    )
  ),
  tags$style(
    HTML(
      "
      * {
        font-family: Roboto, Arial, sans-serif;
        font-weight:400;
      }
      .title {
        text-align: center;
      }
      footer {
        text-align: center;
        padding: 20px 0;
      }
      .btn-default {
        color: white;
        background-color: #f26d80;
        border-color: transparent;
        margin: 0 auto;
        display: block;
      }
      pre#answer,pre#error {
        background: white;
      }
      .tab-content {
        padding-top: 20px;
        min-height: 650px;
      }
      h4{
        font-weight: 400;
      }
      .center {
        font-weight:bold;
        text-align:center;
        margin: 30px auto 0;
      }
      .description {
        padding: 10px 100px 20px;
        text-align:left;
      }
      .description table td, .description table th {
          padding: 10px 20px;
          text-align: center;
      }
      
      .description table th {
          font-weight: bold;
      }
      
      .description table {
          margin: 20px auto;
      }
      .datatables {
        width: 100%!important;
        overflow: overlay;
      }
      .well {
        background-color: #2d3644;
        border: 1px solid #0000003d;
        color:white;
      }
  
      .formatted_calculations {
        background-color: transparent;
        border: none;
        width: 90%;
        position: relative;
        counter-reset: section;
        margin: 0 auto;
        max-height:650px;
        overflow: overlay;
                  
      }
      .calculations_box {
          text-align: center;
          background-color: white;
          border: 1px solid #0000003d;
          padding: 30px 0 0px;
          border-radius: 10px;
          box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.2);
          position: relative;
          overflow: overlay;
      }
      .calculations_box:first-of-type::before {
        counter-set: section;
      }
      .calculations_box:before {
        position: absolute;
        left: 0px;
        top: 0px;
        width: 100%;
        height: 40px;
        background-color: #f26d80;
        color:white;
        border-radius: 9px 9px 0 0;
        font-size: 18px;
        padding: 8px;
        counter-increment: section;
        content: 'Iteration ' counter(section) ': ';
      }
      
      .calculations_box p{
        font-family: 'Courier New', monospace;
        margin: 0 auto;
        width: fit-content;
        text-align: left;
        font-size: 18px;
        font-weight: bold;
      }
      .nav-tabs>li.active>a, .nav-tabs>li.active>a:focus, .nav-tabs>li.active>a:hover {
        color: white;
        cursor: default;
        background-color: #f26d80;
        border: 1px solid #ddd;
        border-bottom-color: transparent;
      }
      ::-webkit-scrollbar {
        
      }
      "
    )
  ),
  titlePanel(div("Numerical Differentiation: Central Difference", class = "title")),

  tabsetPanel(
    tabPanel("Introduction", 
      fluidRow(
        div(
          class = "description",
          h4(class="center",
             "Definition"
          ),
          h4(
            HTML("<strong>Central difference</strong> is used to approximate the derivative of a function at a given point by evaluating the function at nearby points. It's commonly employed in scientific computing, engineering, physics, and other fields where differential equations arise. Ultimately, it provides a <strong>simple and effective way to approximate derivatives numerically</strong>, making it a valuable tool for various applications in numerical analysis and computational science.")
          ),
          h4("The central difference approximation for the first derivative of a function f(x) at a point x is given by:"
          ),
          h4(class="center",
             withMathJax(h4(class = "center", "$$f'(x) = \\frac{f(x+h) - f(x-h)}{2h}$$"))
          ),
          h4("For every iteration:"),
          h4(class="center",
          withMathJax(h4(class = "center", "$$h = \\frac{h}{2^i}$$")
          )),
          
          h4(HTML("Where <i>i</i> is the no. of iterations")),
          
          h4(class="center",
             "Applications"
          ),
          h4(
            HTML("<strong>Solving Partial Differential Equations. (PDEs)</strong> 
                 Central difference schemes are a numerical method used to solve PDEs governing physical phenomena.
                 <br><br><strong>Signal Processing. </strong> 
                 Central differences help in analyzing signals and extracting essential features.
                 <br><br><strong>Computational Fluid Dynamics. </strong>
                 Central differences play a pivotal role in discretizing the Navier-Stokes equations, allowing for numerical simulations of fluid behavior and related phenomena.
                 ")
          ),
          h4(class="center",
            "R Function Formatting Rules"
          ),
          HTML("<table border='1'>
                  <tr>
                    <th>Instead of:</th>
                    <th>Do</th>
                  </tr>
                  <tr>
                    <td>2x</td>
                    <td>2*x</td>
                  </tr>
                  <tr>
                    <td>2(-x)</td>
                    <td>2*(-x)</td>
                  </tr>
                  <tr>
                    <td>sin*x</td>
                    <td>sin(x)</td>
                  </tr>
                  <tr>
                    <td>sin^2(x)</td>
                    <td>sin(x)^2</td>
                  </tr>
                  <tr>
                    <td>logx</td>
                    <td>log10(x)</td>
                  </tr>
                  <tr>
                    <td>logn(x)</td>
                    <td>log(x)</td>
                  </tr>
                  <tr>
                    <td>1*expx</td>
                    <td>1*exp(x)</td>
                  </tr>
                  <tr>
                    <td>absx</td>
                    <td>abs(x)</td>
                  </tr>
              </table>"
          ),
        ),
        align = "center"
      ),
    ),
  # Inputs and results side by side
    tabPanel("Calculate", 
      fluidRow(
        column(
          width = 3,
          sidebarPanel(
            width = 12,
            h4("Calculator"),
            textInput("func", "Function f in R format:", value = "2 * x * cos(2 * x)"),
            numericInput("x_value", "X:", value = 2),
            numericInput("delta_x", "Step size h or Δx:", value = 0.05, step = 0.01),
            numericInput("num_iter", "No. of iterations:", value = 15, min = 1),
          ),
          column(
            width = 12,
            mainPanel(
              width = 12,
              h4(
                style = "font-size: 1.2em;",
                HTML("The approximate derivative of the function at the specified point:")
              ),
              verbatimTextOutput("answer"),
              h4(
                style = "font-size: 1.2em;",
                HTML("Approximation Error:")
              ),
              verbatimTextOutput("error") #what to put here
            )
          )
        ),
        column(
          width = 8,
          tabsetPanel(
            tabPanel("Plot",
                     plotOutput("plot", height = "500px"),
                     h4(
                       style = "font-size: 1.2em; text-align: justify;",
                       HTML(
                         paste(
                           "<div style='max-width: 900px; margin-left: 50px; margin-right: 50px;'>",
                           "This plot displays the <strong>original function</strong>, represented by the solid curve.",
                           "The red dashed line indicates the <strong>tangent line at the specified point</strong>, providing an approximation of the derivative of the function at that point.",
                           "The slope of this tangent line represents the approximate value of the derivative, showing how the function changes with respect to <strong>x</strong>.",
                           "</div>"
                         )
                       )
                     )
            ),
            tabPanel("Steps",
                     uiOutput("formatted_calculations")
            ),
            tabPanel("Table",
                     dataTableOutput("table"),
                     h4(
                       style = "font-size: 1.1em; text-align: justify;",
                       HTML(
                         paste(
                           "<div style='max-width: 900px; margin-left: 20px; margin-right: 20px; margin-top: 30px;'>",
                           "This table showcases the <strong>Central Difference</strong> method for estimating the derivative of a function at a specific point. 
                           Each iteration of our calculation process is displayed in a row. 
                           The <strong>Step Size</strong> column indicates the size of the step taken to approximate the derivative, 
                           while the <strong>Approximated Derivative</strong> column shows our estimation. 
                           The <strong>Approximation Error</strong> measures the difference between our approximation and the actual derivative. 
                           We can refine our approximation for greater accuracy by adjusting parameters such as the function, initial point, step size, and number of iterations. 
                           Generally, decreasing the step size or increasing the number of iterations helps to reduce the approximation error, leading to a more precise estimate.",
                           "</div>"
                         )
                       )
                     )
            )
          )
        )
      )
    )
  ),
  tags$footer(
    strong("Developed by Group 4 (Borromeo, Garcia, Samson)"),
    br(),
    "© 2024", " All rights reserved."
  ),
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive calculation
  df_results <- reactive({
    validate(
      need(input$func != "","Function input is empty"),
      need(input$x_value != "","X input is empty"),
      need(input$delta_x != "","h or delta_X input is empty"),
      need(input$num_iter != "","No. of iterations is empty")
    )
    f <- function(x) {
      eval(parse(text = input$func))
    }
    
    x <- input$x_value
    delta_x <- input$delta_x
    num_iter <- input$num_iter
    
    df <- data.frame(matrix(ncol = 6, nrow = num_iter + 1))
    colnames(df) <- c("Iteration", "Step Size h or (Δx)", "f(x+h)", "f(x-h)", "Approximated Derivative (f'(x))", "Approximation Error")
    
    # Perform central difference calculation for each iteration
    for (i in 0:num_iter) {
      delta_x_i <- delta_x / 2^i
      df$Iteration[i + 1] <- i
      df$`Step Size h or (Δx)`[i + 1] <- delta_x_i
      
      f_x_plus_h <- f(x + delta_x_i)
      f_x_minus_h <- f(x - delta_x_i)
      
      # Add f(x+h) and f(x-h) to the data frame
      df$`f(x+h)`[i + 1] <- f_x_plus_h
      df$`f(x-h)`[i + 1] <- f_x_minus_h
      
      # Calculate approximated derivative
      df$`Approximated Derivative (f'(x))`[i + 1] <- (f_x_plus_h - f_x_minus_h) / (2 * delta_x_i)
      
      # Calculate approximation error
      actual_derivative <- grad(f, x)
      df$`Approximation Error`[i + 1] <- abs(df$`Approximated Derivative (f'(x))`[i + 1] - actual_derivative)
    }
    
    df
  })
  
  # Plot the function and the derivative
  output$plot <- renderPlot({
    df <- df_results()
    f <- function(x) {
      eval(parse(text = input$func))
    }
    x <- input$x_value
    curve(f, -10, 10, main = "Plot of the Function and Approximated Derivative", xlab = "x", ylab = "f(x)",
          cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.2)
    
    delta_x <- input$delta_x
    num_iter <- input$num_iter
    delta_x_i <- delta_x / 2^num_iter
    slope <- central_difference(f, x, delta_x_i)
    
    abline(a = f(x), b = slope, col = "red", lty = 2)
    
    #text(0, 0, paste("Slope of line is f'(",x,") =", slope), pos = 3, col = "red",cex = 1.2)
  })
  
  # Display the final answer
  output$answer <- renderPrint({
    df <- df_results()
    df[nrow(df), "Approximated Derivative (f'(x))"]
  })
  
  output$error <- renderPrint({
    df <- df_results()
    df[nrow(df), "Approximation Error"]
  })
  
  output$table <- renderDataTable({
    df <- df_results()
    datatable(df, options = list(pageLength = 10, searching = FALSE, lengthChange = FALSE), rownames = FALSE)
  })
  
  formatted_calculations_text <- reactive({
    x <- input$x_value
    num_iter <- input$num_iter  
    df <- df_results()  
    
    calculations <- character(nrow(df))
    for (i in 0:num_iter) {
      calculations[i + 1] <- paste0(
        "f(x + h)",
          paste0(" = f(", x, " + ", df$`Step Size h or (Δx)`[i + 1], ") ≈ ", df$`f(x+h)`[i + 1]), "\n",
        "f(x - h)",
        paste0(" = f(", x, " - ", df$`Step Size h or (Δx)`[i + 1], ") ≈ ", df$`f(x-h)`[i + 1]), "\n\n",
        "f'(x)",
        paste0(" = (f(x + h) - f(x - h)) / (2 * ", df$`Step Size h or (Δx)`[i + 1], ") ≈ ", df$`Approximated Derivative (f'(x))`[i + 1]),
        "\n"
      )
    }
    calculations
  })
  
  # Display formatted calculations
  output$formatted_calculations <- renderUI({
    ui <- fluidRow(
      tags$pre(class="formatted_calculations",
               lapply(formatted_calculations_text(), function(text) {
                 div(class = "calculations_box", 
                     h4(text) 
                 )
               })
      )
    )
    ui
  })
  
  
  # Calculate results on button click
  observeEvent(input$calculate, {
    df_results()
  })
  }
  
  # Run app
  shinyApp(ui, server)
  
