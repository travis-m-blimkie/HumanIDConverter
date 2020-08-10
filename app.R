library(shinythemes)
library(shinyjs)
library(shiny)
library(DT)
library(tidyverse)

shiny_biomart_table <- readRDS("data/shiny_biomart_table.rds")


shinyApp(
    ui = fluidPage(

        tags$head(tags$link(
            rel = "stylesheet", type = "text/css", href = "css/user.css"
        )),

        title = "Human Gene ID Converter",

        shinyjs::useShinyjs(),

        theme = "readablebootstrap.css",

        titlePanel(h1(
            "Human Gene ID Conveter",
            style = paste0(
                "background-color: #4582ec; color: white; padding: 20px; ",
                "margin: 0px; margin-bottom: 20px"
            )
        )),

        sidebarLayout(

            sidebarPanel = sidebarPanel(

                tags$p("This is some test text."),

                textAreaInput(
                    inputId = "pastedInput",
                    label   = NULL,
                    placeholder = "Your genes here...",
                    height  = 200
                ),

                actionButton(
                    inputId = "search",
                    label   = "Search",
                    icon    = icon("search")
                )
            ),

            mainPanel = mainPanel(
                tags$br(),
                uiOutput("matchedPanel")
            )
        )
    ),


    server = function(input, output) {

        inputGenes <- reactiveVal()


        observeEvent(input$pastedInput, {
            input$pastedInput %>%
                str_split(., pattern = " |\n") %>%
                unlist() %>%
                inputGenes()
        }, ignoreInit = TRUE, ignoreNULL = TRUE)


        matchedGenes <- reactive({
            req(inputGenes())

            shiny_biomart_table %>%
                filter_all(., any_vars(. %in% inputGenes())) %>%
                distinct(hgnc_symbol, .keep_all = TRUE) %>%
                arrange(hgnc_symbol) %>%
                select(
                    "HGNC" = hgnc_symbol,
                    "Ensembl" = ensembl_gene_id,
                    "Entrez" = entrezgene_id,
                    "UniProt" = uniprot_gn_id
                )
        })


        output$matchedTable <- DT::renderDataTable({
            isolate(matchedGenes())
        },
        rownames = FALSE,
        options = list(
            scrollX = "100%",
            scrollY = "500px",
            scrollCollapse = TRUE,
            paging = FALSE
        ))

        observeEvent(input$search, {
            output$matchedPanel <- renderUI({
                tagList(
                    # tags$br(),
                    DT::dataTableOutput("matchedTable")
                )
            })
        })

    }
)
