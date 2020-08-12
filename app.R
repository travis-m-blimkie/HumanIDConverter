
# Load libraries and data -------------------------------------------------

library(shinythemes)
library(bootstraplib)
library(shinyjs)
library(shiny)
library(DT)
library(tidyverse)

shiny_biomart_table <- readRDS("data/shiny_biomart_table.rds")




# Start the app! ----------------------------------------------------------

shinyApp(


    # First the UI section ------------------------------------------------

    ui = fluidPage(

        # Link to some custom CSS
        tags$head(tags$link(
            rel = "stylesheet", type = "text/css", href = "css/user.css"
        )),

        title = "Human ID Converter",

        shinyjs::useShinyjs(),

        # Select the Bootswatch3 "Readable": https://bootswatch.com/3/readable/
        theme = "readablebootstrap.css",

        # Customize the header for the app, with a theme-appropriate colour
        titlePanel(h1(
            "Human ID Converter",
            style = paste0(
                "background-color: #4582ec; color: white;",
                "padding: 20px; margin: 0px; margin-bottom: 20px"
            )
        )),

        sidebarLayout(

            sidebarPanel = sidebarPanel(

                tags$p(
                    "Hello! This app allows you to input human gene IDs or ",
                    "names, and will search for these in a table and return ",
                    "any matches, facilitating mapping between ID types. ",
                    "Currently we support the following names/ID types: ",
                    "HGNC symbols, Ensembl ID, Entrez IDs, and Uniprot IDs."
                ),

                tags$p(
                    "To get started, paste your IDs or names into the field ",
                    "below (one per line), and click the 'Search' button to ",
                    "begin!"
                ),

                tags$br(),

                # User input box
                textAreaInput(
                    inputId = "pastedInput",
                    label   = NULL,
                    placeholder = "Your genes here...",
                    height  = 300
                ),

                # Search button, which is a trigger for lots of outputs/buttons
                actionButton(
                    inputId = "search",
                    label   = "Search",
                    icon    = icon("search")
                ),

                tags$br(),

                # Download button for matching genes
                uiOutput(outputId = "matchedBtn")
            ),

            mainPanel = mainPanel(

                tags$br(),

                # Output table for matching genes
                uiOutput("matchedPanel")
            )
        )
    ),





    # Now the server section ----------------------------------------------

    server = function(input, output) {

        inputGenes <- reactiveVal()


        # Take the user's input and clean it up
        observeEvent(input$pastedInput, {
            input$pastedInput %>%
                str_split(., pattern = " |\n") %>%
                unlist() %>%
                inputGenes()
        }, ignoreInit = TRUE, ignoreNULL = TRUE)


        # Now look for the user's genes in our table
        matchedGenes <- reactive({
            req(inputGenes())

            shiny_biomart_table %>%
                filter_all(., any_vars(. %in% inputGenes())) %>%
                distinct(hgnc_symbol, .keep_all = TRUE) %>%
                arrange(hgnc_symbol) %>%
                select(
                    "HGNC"    = hgnc_symbol,
                    "Ensembl" = ensembl_gene_id,
                    "Entrez"  = entrezgene_id,
                    "UniProt" = uniprot_gn_id
                )
        })



        # Output tables ---------------------------------------------------

        # First for the matching genes
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
                DT::dataTableOutput("matchedTable")
            })
        })



        # Download buttons ------------------------------------------------

        # First for the matched genes
        output$matchedDl <- downloadHandler(
            filename = "matching_genes.csv",
            content = function(file) {
                readr::write_csv(matchedGenes(), path = file)
            }
        )

        # Render the download button, after the search button is clicked and
        # only when there are some matching genes to return.
        observeEvent(input$search, {
            output$matchedBtn <- renderUI({
                isolate(matchedGenes())

                if (nrow(matchedGenes()) != 0) {
                    tagList(
                        tags$hr(),
                        downloadButton(
                            outputId = "matchedDl",
                            label    = "Matched Genes"
                        ),
                        tags$br()
                    )
                }
            })
        })
    }
)
