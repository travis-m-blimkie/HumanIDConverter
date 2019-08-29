
# Load libraries and global variables -------------------------------------

library(shinythemes)
library(shinyjs)
library(shiny)
library(tidyverse)

shiny_biomart_table <- readRDS("data/shiny_biomart_table.rds")


# Define UI for data upload app -------------------------------------------

ui <- fluidPage(

    theme = shinytheme("flatly"),
    shinyjs::useShinyjs(),

    # App title
    titlePanel(h1("Human Gene ID Converter")),

    # Sidebar layout with input and output definitions
    sidebarLayout(

        # Sidebar panel for inputs
        sidebarPanel(

            # Input: Select a file
            fileInput("file1", paste0("Choose a text file to upload,",
                                      " with one gene per line:"),
                      multiple = TRUE,
                      accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")),

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
                "<b>Search</b> button to retreive your mappings."
            ))),

            tags$hr(),

            # Button which triggers results to display. "background-color"
            # defines the colour of the button (default #337ab7)
            actionButton(
                inputId = "search",
                label   = "Search",
                icon    = icon("search"),
                style   = "color: #fff; background-color: #18bc9c; border-color: #18bc9c; width: 130px"
            ),

            tags$hr(),

            # Download matched genes
            disabled(downloadButton(
                "Matched_dl",
                "Download Matched Genes",
                style = "width: 260px; background-color: #2c3e50; border-color: #2c3e50"
            )),

            tags$br(),
            tags$br(),

            # Download non-matching genes
            disabled(downloadButton(
                "NoMatch_dl",
                "Download Non-Matching Genes",
                style = "width: 260px; background-color: #2c3e50; border-color: #2c3e50"
            )),

            tags$br(),
            tags$br(),

            # Download "LOC" genes
            disabled(downloadButton(
                "LOC_dl",
                "Download LOC Genes",
                style = "width: 260px; background-color: #2c3e50; border-color: #2c3e50"
            ))

        ),

        # Main panel for displaying results
        mainPanel(

            # Output for matching genes
            h3("Matching Genes:\n"),
            tableOutput("matchedGenes"),

            # Output for non-matching genes
            h3("Non-matching Genes:\n"),
            tableOutput("nonMatchedGenes"),

            # Output for "LOC" genes
            h3("LOC Genes:\n"),
            tableOutput("LOCGenes")

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

    observeEvent(input$file1, {
        showNotification("File successfully uploaded!", type = "message")
    })

    observeEvent(input$search, {

        # Clean input genes ---------------------------------------------------

        # Remove leading or trailing spaces, and remove genes containing a "."
        cleanGenes_1 <- reactive({
            req(inputGenes)

            inputGenes()[, 1] %>%
                as.character() %>%
                str_trim(., ) %>%
                str_subset(., pattern = "\\.", negate = TRUE)
        })


        # Remove genes containing LOC -----------------------------------------

        cleanGenes_2 <- reactive({
            req(cleanGenes_1)

            grep(cleanGenes_1(),
                 pattern = "^LOC[0-9]+$",
                 value = TRUE,
                 invert = TRUE
            )
        })


        # Look for matching entries -------------------------------------------

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


        # Find genes in original list with no matches -------------------------

        nonMatchedGenes <- reactive({
            req(cleanGenes_2, matchedGenes)

            myMatches <- unlist(matchedGenes()) %>% as.character()
            noMatches <- data.frame(Genes = setdiff(cleanGenes_2(), myMatches))

            return(noMatches)
        })


        # Grab genes starting with LOC ----------------------------------------

        LOCGenes <- reactive({
            req(cleanGenes_1)

            LOCs <- grep(cleanGenes_1(),
                         pattern = "^LOC[0-9]+$",
                         value = TRUE)

            data.frame(Genes = LOCs)
        })


        # Display outputs ---------------------------------------------------------

        # Matching Genes
        output$matchedGenes <- renderTable({
            isolate(matchedGenes())
        },
        striped = TRUE
        )

        # Genes without a match
        output$nonMatchedGenes <- renderTable({
            isolate(nonMatchedGenes())
        },
        striped = TRUE)

        # LOC genes which were excluded
        output$LOCGenes <- renderTable({
            isolate(LOCGenes())
        },
        striped = TRUE)


        # Download links ------------------------------------------------------

        enable("Matched_dl")
        output$Matched_dl <- downloadHandler(
            filename = "matching_genes.csv",
            content = function(file) {
                write.csv(matchedGenes(), file, row.names = FALSE, quote = FALSE)
            }
        )

        enable("NoMatch_dl")
        output$NoMatch_dl <- downloadHandler(
            filename = "non-matching_genes.csv",
            content = function(file) {
                write.csv(nonMatchedGenes(), file, row.names = FALSE, quote = FALSE)
            }
        )

        enable("LOC_dl")
        output$LOC_dl <- downloadHandler(
            filename = "loc_genes.csv",
            content = function(file) {
                write.csv(LOCGenes(), file, row.names = FALSE, quote = FALSE)
            }
        )

    })


}




# Run the app -------------------------------------------------------------

shinyApp(ui, server)
