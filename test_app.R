

# Minimal version of the app to try and work out how to properly set up reactive
# values and observations to get desired behaviour.


library(shiny)
library(DT)
library(tidyverse)


biomart_table <- readRDS("data/biomart_table.Rds")


ui <- fluidPage(
  sidebarLayout(
    sidebarPanel = sidebarPanel(


      textAreaInput(
        inputId     = "pastedInput",
        label       = NULL,
        placeholder = "Your genes here...",
        height = 200
      ),


      actionButton(
        inputId = "search",
        label   = "Search Genes",
        icon    = icon("search")
      )
    ),

    mainPanel = mainPanel(
      uiOutput("matchedPanel"),
    )
  )
)

server <- function(input, output, session) {


  inputGenes <- reactiveVal()

  observeEvent(input$search, {
    input$pastedInput %>%
      read_lines() %>%
      inputGenes()
  }, ignoreInit = TRUE, ignoreNULL = TRUE)


  matchedGenes <- reactive({
    req(inputGenes)

    biomart_table %>%
      filter_all(., any_vars(. %in% inputGenes())) %>%
      distinct(HGNC, .keep_all = TRUE) %>%
      arrange(HGNC)
  })


  output$matchedTable <- DT::renderDataTable(
    matchedGenes(),
    rownames  = FALSE,
    escape    = FALSE,
    selection = "none",
    options   = list(
      scrollX = "100%",
      scrollY = "250px",
      scrollCollapse = TRUE,
      paging  = FALSE
    )
  )

  output$matchedPanel <- renderUI({
    req(matchedGenes())

    if (nrow(matchedGenes()) != 0) {
      return(
        tagList(
          tags$h3("We found matches for the following genes:"),
          DT::dataTableOutput("matchedTable")
        )
      )
    } else {
      return(NULL)
    }
  })


}

shinyApp(ui, server)
