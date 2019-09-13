
# Load libraries and global variables -------------------------------------

library(shinythemes)
library(shinyjs)
library(shiny)
library(DT)
library(tidyverse)

shiny_biomart_table <- readRDS("data/shiny_biomart_table.rds")

# Colours that match with the Yeti theme:
#     --blue: #008cba;
#     --purple: #6f42c1;
#     --pink: #e83e8c;
#     --red: #F04124;
#     --orange: #fd7e14;
#     --yellow: #E99002;
#     --green: #43ac6a;
#     --teal: #20c997;
#     --cyan: #5bc0de;


# Define UI for data upload app -------------------------------------------

ui <- fluidPage(

    # Set title that shows in the browser
    title = "Human ID Converter",

    # Set theme
    theme = shinytheme("yeti"),

    # Enable shinyjs usage
    shinyjs::useShinyjs(),

    # App title
    titlePanel(h1("Human Gene ID Converter")),

    # Sidebar layout with input and output definitions
    sidebarLayout(

        # Sidebar panel for inputs
        sidebarPanel(

            # Input: Select a file
            fileInput(
                placeholder = " Please select a file",
                "file1",
                paste0("Choose a text file to upload,",
                       " with one gene per line:"),
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")
            ),

            tags$hr(),

            # Select separator for input list (not really needed if one gene per
            # line)
            radioButtons("sep", "Separator",
                         choices = c(Comma = ",",
                                     Semicolon = ";",
                                     Tab = "\t"),
                         selected = ","),


            # Input: Checkbox if file has header
            checkboxInput("header", "File contains header", FALSE),

            tags$hr(),

            tags$p(div(HTML(
                "Once you're genes have been successfully uploaded, press the",
                "<b>Search</b> button to retreive your mappings and download ",
                "the results."
            ))),

            tags$hr(),

            # Search button which triggers results to display
            actionButton(
                inputId = "search",
                label   = "Search",
                icon    = icon("search"),
                style   = "width: 130px; background-color: #008cba; border-color: #008cba; color: #ffffff"
            ),

            tags$hr(),

            # All the buttons below (Matched, NonMatched, LOCGenes) are created
            # with uiOutput() and renderUI() for greater control.

            # Download matched genes
            uiOutput("Matched_btn"),


            # Download non-matching genes
            uiOutput("NonMatched_btn"),


            # Download "LOC" genes
            uiOutput("LOCGenes_btn")

        ),

        # Main panel for displaying results
        mainPanel(

            # Output for matching genes, just rendered plainly since this will
            # (presumably) always contain something
            h3("Your matching genes will be displayed below:"),
            DT::dataTableOutput("matchedGenes"),

            # Output for non-matching genes, using uiOutput() and renderUI() so
            # it only displays when populated
            uiOutput("NonMatchingGenes_tbl"),

            # Output for "LOC" genes, same as NonMatchingGenes_tbl
            uiOutput("LOCGenes_tbl"),

            tags$hr()

        )

    )
)




# Define server logic to read selected file -----------------------------------

server <- function(input, output) {


    # Read input genes from user ------------------------------------------

    inputGenes <- reactive({
        req(input$file1)

        read.csv(input$file1$datapath,
                 header = input$header,
                 sep = input$sep)
    })


    # Notification for successful upload ----------------------------------

    # observeEvent(input$file1, {
    #     showNotification("File successfully uploaded!", type = "message")
    # })

    # Start of main server code -------------------------------------------
    observeEvent(input$search, {

        # Clean input genes -----------------------------------------------

        # Remove leading or trailing spaces, and remove genes containing a "."
        cleanGenes_1 <- reactive({
            req(inputGenes)

            inputGenes()[, 1] %>%
                as.character() %>%
                str_trim(., ) %>%
                str_subset(., pattern = "\\.", negate = TRUE)
        })


        # Remove genes containing LOC -------------------------------------

        cleanGenes_2 <- reactive({
            req(cleanGenes_1)

            grep(cleanGenes_1(),
                 pattern = "^LOC[0-9]+$",
                 value = TRUE,
                 invert = TRUE
            )
        })


        # Look for matching entries ---------------------------------------

        matchedGenes <- reactive({
            req(cleanGenes_2)

            shiny_biomart_table %>%
                filter_all(., any_vars(. %in% cleanGenes_2())) %>%
                distinct(hgnc_symbol, .keep_all = TRUE) %>%
                arrange(hgnc_symbol) %>%
                rename(
                    "HGNC" = hgnc_symbol,
                    "UniProt" = uniprot_gn_id,
                    "Entrez" = entrezgene_id,
                    "Ensembl" = ensembl_gene_id
                )
        })


        # Find genes in original list with no matches ---------------------

        nonMatchedGenes <- reactive({
            req(cleanGenes_2, matchedGenes)

            myMatches <- unlist(matchedGenes()) %>% as.character()
            noMatches <- data.frame(Genes = setdiff(cleanGenes_2(), myMatches))

            return(noMatches)
        })


        # Grab genes starting with LOC ------------------------------------

        LOCGenes <- reactive({
            req(cleanGenes_1)

            LOCs <- grep(cleanGenes_1(),
                         pattern = "^LOC[0-9]+$",
                         value = TRUE)

            data.frame(Genes = LOCs)
        })


        # Display outputs -------------------------------------------------

        # Matching Genes, rendered simply
        output$matchedGenes <- DT::renderDataTable({
            isolate(matchedGenes())
        }, options = list(searching = FALSE,
                           scrollX = "100%",
                           scrollY = "400px",
                           scrollCollapse = TRUE,
                           paging = FALSE),
        rownames = FALSE
        )

        # Genes without a match, making use of renderUI()
        output$nonMatchedGenes_DT <- DT::renderDataTable({
            isolate(nonMatchedGenes())
        }, options = list(searching = FALSE,
                          scrollX = "100%",
                          scrollY = "250px",
                          scrollCollapse = TRUE,
                          paging = FALSE),
        rownames = FALSE
        )

        # Rendering the non-matched genes if present
        output$NonMatchingGenes_tbl <- renderUI({
            isolate(nonMatchedGenes())

            if (nrow(nonMatchedGenes()) == 0) {
                return(NULL)
            } else {
                return(tagList(
                    tags$hr(),

                    tags$h3("These genes had no matches in our database:"),
                    DT::dataTableOutput("nonMatchedGenes_DT")
                ))
            }
        })

        # LOC genes which were excluded
        output$LOCGenes_DT <- DT::renderDataTable({
            isolate(LOCGenes())
        }, options = list(searching = FALSE,
                          scrollX = "100%",
                          scrollY = "250px",
                          scrollCollapse = TRUE,
                          paging = FALSE),
        rownames = FALSE
        )

        # Rendering the LOC genes if present
        output$LOCGenes_tbl <- renderUI({
            isolate(LOCGenes())

            if (nrow(LOCGenes()) == 0) {
                return(NULL)
            } else {
                return(tagList(
                    tags$hr(),

                    tags$h3("The following LOC genes were detected:"),
                    DT::dataTableOutput("LOCGenes_DT")
                ))
            }
        })


        # Download buttons ------------------------------------------------

        # First for the matched genes -------
        output$Matched_dl <- downloadHandler(
            filename = "matching_genes.csv",
            content = function(file) {
                write.csv(matchedGenes(), file, row.names = FALSE, quote = FALSE)
            }
        )

        output$Matched_btn <- renderUI({
            isolate(matchedGenes())

            if(nrow(matchedGenes()) != 0) {
                tagList(
                    downloadButton(
                        "Matched_dl",
                        "Matched Genes",
                        style = "width: 200px"
                    ),
                    tags$br(),
                    tags$br()
                )
            }
        })


        # Next for non-matching genes -------
        output$NoMatch_dl <- downloadHandler(
            filename = "non-matching_genes.csv",
            content = function(file) {
                write.csv(nonMatchedGenes(), file, row.names = FALSE, quote = FALSE)
            }
        )

        output$NonMatched_btn <- renderUI({
            isolate(nonMatchedGenes())

            if (nrow(nonMatchedGenes()) != 0) {
                tagList(
                    downloadButton(
                        "NoMatch_dl",
                        "Non-Matching Genes",
                        style = "width: 200px"
                    ),
                    tags$br(),
                    tags$br()
                )
            }
        })


        # Finally for LOC genes ---------
        output$LOC_dl <- downloadHandler(
            filename = "loc_genes.csv",
            content = function(file) {
                write.csv(LOCGenes(), file, row.names = FALSE, quote = FALSE)
            }
        )

        output$LOCGenes_btn <- renderUI({
            isolate(LOCGenes())

            if (nrow(LOCGenes()) != 0) {
                tagList(
                    downloadButton(
                        "LOC_dl",
                        "LOC Genes",
                        style = "width: 200px"
                    )
                )
            }
        })

    })


}




# Run the app -------------------------------------------------------------

shinyApp(ui, server)
