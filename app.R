
# Load libraries and data -------------------------------------------------

library(shiny)
library(DT)
library(tidyverse)

# These two lines are needed if publishing the app the RStudio's shiny service,
# so the packages from BioConductor can be found and installed
# library(BiocManager)
# options(repos = BiocManager::repositories())

biomart_table <- readRDS("data/biomart_table.Rds")
example_data  <- read_lines("data/shiny_app_test_data.txt")




# First define the UI section ---------------------------------------------

ui <- fluidPage(

    title = "Human ID Converter",

    # Link to some custom CSS tweaks
    tags$head(tags$link(
        rel = "stylesheet", type = "text/css", href = "css/user.css"
    )),

    # Select the Bootswatch3 theme "Readable": https://bootswatch.com/3/readable
    theme = "css/readablebootstrap.css",

    # Header/title of the app, which has some custom tweaks applied in
    # "www/css/user.css"
    tags$h1("Human ID Converter"),

    sidebarLayout(

        sidebarPanel = tags$div(

            class = "col-sm-4 manual-sidebar",

            tags$form(
                class = "well",

                tags$h3("Hello!"),

                tags$p(HTML(
                    "Welcome to <b>Human ID Converter</b>, an app ",
                    "designed to facilitate mapping between different ",
                    "human gene identifiers. It functions by searching a ",
                    "table with your input genes and returning any ",
                    "matches. Each input can contain a mix of any of the ",
                    "supported ID types: HGNC, Ensembl, & Entrez."
                )),

                tags$p(
                    "The data used for the mapping comes from Ensembl's ",
                    tags$a(
                        "BioMart",
                        href = "http://ensemblgenomes.org/info/access/biomart",
                        .noWS = c("before", "after")
                    ),
                    ". If you run into any trouble, please ",
                    "open an issue at the ",
                    tags$a(
                        "Github page",
                        href = "https://github.com/travis-m-blimkie/HumanIDConverter",
                        .noWS = c("before", "after"),
                    ),
                    "."
                ),

                tags$p(
                    "To get started, paste your genes into the field ",
                    "below (one per line), and click the 'Search Genes' ",
                    "button to see your results."
                ),

                tags$br(),

                # Field for user to input their genes
                textAreaInput(
                    inputId     = "pastedInput",
                    label       = NULL,
                    placeholder = "Your genes here...",
                    height      = 155
                ),

                # Link to load example data, primarily to making testing
                # easier
                actionLink(
                    inputId = "tryExample",
                    label   = "Load Example Data",
                    style   = "font-size: 110%"

                ),

                # Search button, which is a trigger for lots of
                # outputs/buttons
                actionButton(
                    class   = "btn-primary",
                    style   = "float: right; padding-bottom: 10px",
                    inputId = "search",
                    label   = "Search Genes",
                    icon    = icon("search")
                ),

                tags$br(),
                tags$br(),

                # Download button and some text for matching genes
                uiOutput("matchedBtn"),

                # Download button and some text for non-matching genes
                uiOutput("nonMatchedBtn")

            ) # Closes the well/form

        ), # Closes sidebar div/well

        mainPanel = tags$div(

            class = "col-sm-8",

            tags$br(),

            # Output table of matching genes
            uiOutput("matchedPanel"),

            # Output table of non-matching genes
            uiOutput("nonMatchedPanel")

        ) # Closes mainPanel() div

    ) # Closes the sidebarLayout()

) # Closes the fluidPage() UI call





# Now the server section --------------------------------------------------


server <- function(input, output) {

    inputGenes <- reactiveVal()


    # Load in example data when linked clicked, and provide a notification. Note
    # the "message" notification type has been modified; see "www/css/user.css"
    # for details.
    observeEvent(input$tryExample, {
        showNotification(
            id = "exampleSuccess",
            duration = 5,
            closeButton = TRUE,
            type = "message",
            ui = paste0(
                "Example data successfully loaded! Click the 'Search Genes' ",
                "button to continue."
            )
        )
    })

    observeEvent({
        input$tryExample
        input$search
    }, {
        inputGenes(example_data)
    })


    # Take the user's input and clean it up, matching a space or new line
    # character to separate out the genes into a character vector
    observeEvent(input$search, {
        input$pastedInput %>%
            str_split(., pattern = " |\n") %>%
            unlist() %>%
            inputGenes()
    }, ignoreInit = TRUE, ignoreNULL = TRUE)


    # Now look for the user's genes in the biomaRt table. By using
    # `filter_all()` and `any_vars()` we can search every column of the
    # biomaRt table for each input gene, meaning the user can input a mixture
    # of different ID types and we don't need specific code for each.
    hyperlinkConstructor <- reactive({
        req(inputGenes())

        biomart_table %>%
            filter_all(., any_vars(. %in% inputGenes())) %>%
            distinct(HGNC, .keep_all = TRUE) %>%
            arrange(HGNC)
    })


    matchedGenes <- reactive({
        req(hyperlinkConstructor())
        hyperlinkConstructor() %>% select(-HGNC_ID)
    })



    # Create an alternate table to display to the user, in which all the genes
    # are links to the respective info page for that gene.
    # Note for HGNC we use the corresponding HGNC ID to link to the proper page.
    # This column gets dropped at the end, since it's just used for constructing
    # the link.
    hyperlinkTable <- reactive({
        req(hyperlinkConstructor())

        hyperlinkConstructor() %>%
            rowwise() %>%
            mutate(
                HGNC = paste0(
                    "<a target='_blank' href='",
                    "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:",
                    HGNC_ID, "'>", HGNC, "</a>"
                ),
                Ensembl = paste0(
                    "<a target='_blank' href='",
                    "http://www.ensembl.org/id/",
                    Ensembl, "'>", Ensembl, "</a>"
                ),
                Entrez = paste0(
                    "<a target='_blank' href='",
                    "https://www.ncbi.nlm.nih.gov/gene/",
                    Entrez, "'>", Entrez, "</a>"
                ),
                UniProt = paste0(
                    "<a target='_blank' href='",
                    "https://www.uniprot.org/uniprot/",
                    UniProt, "'>", UniProt, "</a>"
                )
            ) %>%
            ungroup() %>%
            select(-HGNC_ID)
    })

    # Get the genes that didn't have matches to inform the user
    nonMatchedGenes <- reactive({
        req(matchedGenes())

        myMatches <- unlist(matchedGenes()) %>% as.character()
        noMatches <-
            tibble("Input Genes" = setdiff(inputGenes(), myMatches)) %>%
            arrange(`Input Genes`)

        return(noMatches)
    })

    # Create a character vector version of non-matching genes for download
    # purposes (we still want the above tibble version with links for displaying
    # to the user).
    nonMatchedGenes_chr <- reactive({
        req(nonMatchedGenes())
        nonMatchedGenes() %>% pull(`Input Genes`)
    })



    # Output tables -------------------------------------------------------

    # First for the matching genes
    output$matchedTable <- DT::renderDataTable(
        hyperlinkTable(),
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
        req(hyperlinkTable())

        if (nrow(hyperlinkTable()) != 0) {
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





    # Now for non-matching genes
    output$nonMatchedTable <- DT::renderDataTable(
        nonMatchedGenes(),
        rownames  = FALSE,
        selection = "none",
        options   = list(
            scrollX = "100%",
            scrollY = "250px",
            scrollCollapse = TRUE,
            paging  = FALSE
        )
    )

    output$nonMatchedPanel <- renderUI({
        req(nonMatchedGenes())

        if (nrow(nonMatchedGenes()) != 0) {
            return(tagList(
                tags$hr(),
                tags$br(),
                tags$h3("The following genes did not have a match:"),
                DT::dataTableOutput("nonMatchedTable"),
                tags$br()
            ))

        } else {
            return(NULL)
        }
    })






    # Download buttons ----------------------------------------------------

    # First for the matched genes
    output$matchedDl <- downloadHandler(
        filename = "matching_genes.csv",
        content  = function(f) {
            readr::write_csv(matchedGenes(), file = f)
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
                        "the table."
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
        filename = "non_matching_genes.txt",
        content  = function(f) {
            readr::write_lines(nonMatchedGenes_chr(), file = f)
        }
    )

    observeEvent(input$search, {
        output$nonMatchedBtn <- renderUI({
            isolate(nonMatchedGenes())

            if (nrow(nonMatchedGenes()) != 0) {
                tagList(
                    tags$hr(),
                    tags$p(
                        "We were unable to find matches for the genes ",
                        "shown in the bottom table on the right. Click ",
                        "the button below to download them."
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




# Call the app to start ---------------------------------------------------

shinyApp(ui = ui, server = server)
