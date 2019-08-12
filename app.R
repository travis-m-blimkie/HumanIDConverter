
# Load libraries and global variables -------------------------------------

library(tidyverse)
library(shinythemes)
library(shiny)

shiny_biomart_table <- readRDS("data/shiny_biomart_table.rds")


# Define UI for data upload app -------------------------------------------

ui <- fluidPage(

    theme = shinytheme("flatly"),

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

            # Download matched genes
            downloadButton("Matched_dl", "Download Matched Genes"),
            tags$br(),
            tags$br(),

            # Download non-matching genes
            downloadButton("NoMatch_dl", "Download Non-Matching Genes"),
            tags$br(),
            tags$br(),

            # Download "LOC" genes
            downloadButton("LOC_dl", "Download LOC Genes")

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
    output$matchedGenes <- renderTable(
        {matchedGenes()},
        striped = TRUE
        )

    # Genes without a match
    output$nonMatchedGenes <- renderTable({
        nonMatchedGenes()
    },
    striped = TRUE)

    # LOC genes which were excluded
    output$LOCGenes <- renderTable({
        LOCGenes()
    },
    striped = TRUE)


    # Download links ------------------------------------------------------
    output$matchedGenes_dl <- downloadHandler(
        filename = "matching_genes.csv",
        content = function(file) {
            write.csv(matchedGenes(), file, row.names = FALSE, quote = FALSE)
        }
    )

    output$nonMatchedGenes_dl <- downloadHandler(
        filename = "non-matching_genes.csv",
        content = function(file) {
            write.csv(nonMatchedGenes(), file, row.names = FALSE, quote = FALSE)
        }
    )

    output$LOCGenes_dl <- downloadHandler(
        filename = "loc_genes.csv",
        content = function(file) {
            write.csv(LOCGenes(), file, row.names = FALSE, quote = FALSE)
        }
    )



}




# Run the app -------------------------------------------------------------

shinyApp(ui, server)
