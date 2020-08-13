
# Load libraries and data -------------------------------------------------

library(shiny)
library(DT)
library(tidyverse)

shiny_biomart_table <- readRDS("data/shiny_biomart_table.rds")




# Start the app! ----------------------------------------------------------

shinyApp(


    # First the UI section ------------------------------------------------

    ui = fluidPage(

        # Link to some custom CSS tweaks
        tags$head(tags$link(
            rel = "stylesheet", type = "text/css", href = "user.css"
        )),

        title = "Human ID Converter",

        # Select the Bootswatch3 "Readable": https://bootswatch.com/3/readable/
        theme = "readablebootstrap.css",

        tags$h1("Human ID Converter"),

        sidebarLayout(

            sidebarPanel = sidebarPanel(

                tags$p(
                    "Hello! This app allows you to input human gene IDs or ",
                    "names, and will search for these in a table and return ",
                    "any matches, facilitating mapping between ID types. ",
                    "Currently we support the following names/ID types: ",
                    "HGNC symbols, Ensembl ID, Entrez IDs, and Uniprot IDs. ",
                    "Data for this app comes from Ensembl's ",
                    tags$a(
                        href = "http://ensemblgenomes.org/info/access/biomart",
                        "BioMart.",
                        .noWS = c("before", "after")
                    )
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
                    height  = 200
                ),

                # Search button, which is a trigger for lots of outputs/buttons
                actionButton(
                    class = "btn btn-primary",
                    style = "float: right; padding-bottom: 10px",
                    inputId = "search",
                    label   = "Search",
                    icon    = icon("search")
                ),

                tags$br(),
                tags$br(),

                # Download button for matching genes
                uiOutput("matchedBtn"),

                uiOutput("nonMatchedBtn")
            ),

            mainPanel = mainPanel(

                tags$br(),

                # Output table for matching genes
                uiOutput("matchedPanel"),

                uiOutput("nonMatchedPanel")
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

        # Get the genes that didn't have matches
        nonMatchedGenes <- reactive({
            req(matchedGenes())

            myMatches <- unlist(matchedGenes()) %>% as.character()
            noMatches <- tibble("Input Genes" = setdiff(inputGenes(), myMatches))

            return(noMatches)
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
        }, ignoreNULL = TRUE, ignoreInit = TRUE)


        # Now for non-matching genes
        output$nonMatchedTable <- DT::renderDataTable({
            isolate(nonMatchedGenes())
        },
        rownames = FALSE,
        options = list(
            scrollX = "100%",
            scrollY = "500px",
            scrollCollapse = TRUE,
            paging = FALSE
        ))

        observeEvent(input$search, {
            output$nonMatchedPanel <- renderUI({
                isolate(nonMatchedGenes())

                if (nrow(nonMatchedGenes()) == 0) {
                    return(NULL)
                } else {
                    return(tagList(
                        tags$hr(),
                        tags$br(),
                        tags$h3("The following genes did not have a match:"),
                        DT::dataTableOutput("nonMatchedTable"),
                        tags$br()
                    ))
                }
            })
        }, ignoreNULL = TRUE, ignoreInit = TRUE)



        # Download buttons ------------------------------------------------

        # First for the matched genes
        output$matchedDl <- downloadHandler(
            filename = "matching_genes.csv",
            content = function(file) {
                readr::write_csv(matchedGenes(), path = file)
            }
        )

        observeEvent(input$search, {
            output$matchedBtn <- renderUI({
                isolate(matchedGenes())

                if (nrow(matchedGenes()) != 0) {
                    tagList(
                        tags$br(),
                        tags$hr(),
                        tags$p(
                            "Success! We found matches for some of your input ",
                            "genes. Check the table on the right to see which ",
                            "genes we were able to identify. Click the button ",
                            "below to download a CSV table of your result:"
                        ),
                        downloadButton(
                            class    = "btn btn-success",
                            style    = "float: right; padding-bottom: 10px",
                            outputId = "matchedDl",
                            label    = "Matched Genes"
                        ),
                        tags$br(),
                        tags$br()
                    )
                }
            })
        }, ignoreNULL = TRUE, ignoreInit = TRUE)

        output$nonMatchedDl <- downloadHandler(
            filename = "non_matching_genes.csv",
            content = function(file) {
                readr::write_csv(nonMatchedGenes(), path = file)
            }
        )

        # Now for the non-matching genes
        observeEvent(input$search, {
            output$nonMatchedBtn <- renderUI({
                isolate(nonMatchedGenes())

                if (nrow(nonMatchedGenes()) != 0) {
                    tagList(
                        tags$br(),
                        tags$hr(),
                        tags$p(
                            "It seems like we were unable to find matches for ",
                            "some of the genes you submitted. See the bottom ",
                            "table on the right to check which genes were not ",
                            "matched in our database. You may also download ",
                            "these genes in a text file using the button below:"
                        ),
                        downloadButton(
                            class    = "btn btn-warning",
                            style    = "float: right; padding-bottom: 10px",
                            outputId = "nonMatchedDl",
                            label    = "Non-Matching Genes"
                        ),
                        tags$br(),
                        tags$br()
                    )
                }
            })
        }, ignoreNULL = TRUE, ignoreInit = TRUE)
    }
)
