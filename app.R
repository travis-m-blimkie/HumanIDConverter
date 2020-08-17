
# Load libraries and data -------------------------------------------------

suppressPackageStartupMessages({
    library(shiny)
    library(DT)
    library(tidyverse)
})

shiny_biomart_table <- readRDS("data/shiny_biomart_table_20200814.Rds")

example_data <-
    read_csv(
        "example_data/shiny_app_test_data.csv",
        col_names = "gene",
        col_types = cols()
    ) %>%
    pull(1)




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

        # Header/title of the app, which has some custom tweaks applied in
        # "user.css"
        tags$h1("Human ID Converter"),

        sidebarLayout(

            sidebarPanel = sidebarPanel(

                tags$h3("Hello!"),

                tags$p(HTML(
                    "Welcome to <b>Human ID Converter</b>, an app designed to ",
                    "facilitate mapping between different human gene ",
                    "identifiers. It functions by searching a table ",
                    "with your input genes and returning any matches. ",
                    "Each input can contain a mix of any of the supported ID ",
                    "types: HGNC, Ensembl, Entrez, and UniProt."
                )),

                tags$p(
                    "The data used for the mapping comes from Ensembl's ",
                    tags$a(
                        href = "http://ensemblgenomes.org/info/access/biomart",
                        "BioMart.",
                        .noWS = "before"
                    ),
                    "If you run into any trouble, please open an issue at the ",
                    tags$a(
                        "Github page.",
                        .noWS = c("before", "after"),
                        href = "https://github.com/travis-m-blimkie/HumanIDConverter"
                    )
                ),

                tags$p(
                    "To get started, paste your genes into the field ",
                    "below (one per line), and click the 'Search' button to ",
                    "see your results."
                ),

                tags$br(),

                # Field for user to input their genes
                textAreaInput(
                    inputId = "pastedInput",
                    label   = NULL,
                    placeholder = "Your genes here...",
                    height  = 175
                ),

                # Link to load example data, primarily to making testing easier
                actionLink(
                    inputId = "tryExample",
                    label   = "Load Example Data",
                    style   = "font-size: 110%"

                ),

                # Search button, which is a trigger for lots of outputs/buttons
                actionButton(
                    class = "btn-primary",
                    style = "float: right; padding-bottom: 10px",
                    inputId = "search",
                    label   = "Search",
                    icon    = icon("search")
                ),

                tags$br(),
                tags$br(),

                # Download button and some text for matching genes
                uiOutput("matchedBtn"),

                # Download button and some text for non-matching genes
                uiOutput("nonMatchedBtn")

            ), # Closes sidebarPanel()

            mainPanel = mainPanel(
                tags$br(),

                # Output table of matching genes
                uiOutput("matchedPanel"),

                # Output table of non-matching genes
                uiOutput("nonMatchedPanel")

            ) # Closes mainPanel()

        ) # Closes sidebarLayout()

    ), # Closes the ui() call





    # Now the server section ----------------------------------------------

    server = function(input, output) {

        inputGenes <- reactiveVal()

        # Load in example data when linked clicked
        observeEvent(input$tryExample, {
            inputGenes(example_data)

            showNotification(
                id = "exampleSuccess",
                ui = "Example data successfully loaded! Click the 'Search' ",
                "button to continue.",
                duration = 5,
                closeButton = TRUE,
                type = "message"
            )
        })


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
            noMatches <- tibble(
                "Input Genes" = setdiff(inputGenes(), myMatches)
            ) %>%
                arrange(1)

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
            scrollY = "250px",
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
            scrollY = "250px",
            scrollCollapse = TRUE,
            paging  = FALSE
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
                        tags$hr(),
                        tags$p(
                            "We successfully matched some of your input ",
                            "genes. Check the table on the right to see the ",
                            "results, and click the button below to download ",
                            "the result."
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




        # Now for the non-matching genes
        output$nonMatchedDl <- downloadHandler(
            filename = "non_matching_genes.csv",
            content = function(file) {
                readr::write_csv(nonMatchedGenes(), path = file)
            }
        )

        observeEvent(input$search, {
            output$nonMatchedBtn <- renderUI({
                isolate(nonMatchedGenes())

                if (nrow(nonMatchedGenes()) != 0) {
                    tagList(
                        tags$hr(),
                        tags$p(
                            "We were unable to find matches for some of your ",
                            "genes, as shown in the bottom table on the right.",
                            " Click the button below to download them."
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

    } # Closes the server() call
)
