library(shiny)
library(reticulate)
library(shinythemes)
library(DT)
library(plotly)

# Configure Python environment
use_python("/usr/bin/python3")
source_python("ivablib/case_report_analyzer.py")

# Initialize analyzer
analyzer <- CaseReportAnalyzer()

# UI Definition
ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("Ivabradine TdP Case Report Analyzer"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("drug_name", "Drug Name", value = "ivabradine"),
      actionButton("analyze_btn", "Analyze", class = "btn-primary"),
      hr(),
      helpText("This tool analyzes PubMed case reports for potential TdP risks.")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Results",
                 h3("Analysis Results"),
                 verbatimTextOutput("analysis_output"),
                 h3("Case Reports"),
                 DTOutput("cases_table")
        ),
        tabPanel("About",
                 h3("About this Tool"),
                 p("This application analyzes PubMed case reports to identify potential Torsades de Pointes (TdP) risks associated with medications."),
                 p("It uses natural language processing to analyze case reports and identify relevant clinical features.")
        )
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  
  # Reactive value for analysis results
  results <- reactiveVal(NULL)
  
  # Analyze button handler
  observeEvent(input$analyze_btn, {
    # Show progress
    withProgress(message = 'Analyzing...', value = 0, {
      
      incProgress(0.3, detail = "Searching PubMed")
      
      # Run analysis
      tryCatch({
        result <- analyzer$analyze_drug(input$drug_name)
        results(result)
        
        incProgress(0.7, detail = "Processing results")
        
      }, error = function(e) {
        showNotification(
          paste("Error:", e$message),
          type = "error"
        )
      })
    })
  })
  
  # Output for analysis results
  output$analysis_output <- renderPrint({
    req(results())
    results()
  })
  
  # Output for cases table
  output$cases_table <- renderDT({
    req(results())
    
    # Extract case reports from results and create data frame
    cases <- results()$case_reports
    if (length(cases) > 0) {
      cases_df <- data.frame(
        Title = sapply(cases, function(x) x$title),
        Authors = sapply(cases, function(x) paste(x$authors, collapse = ", ")),
        Year = sapply(cases, function(x) x$year),
        PMID = sapply(cases, function(x) x$pmid)
      )
      
      datatable(cases_df,
                options = list(pageLength = 5),
                rownames = FALSE)
    }
  })
}

# Run the app
shinyApp(ui = ui, server = server)
