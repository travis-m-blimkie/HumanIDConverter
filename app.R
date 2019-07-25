
# Load libraries, global variables ----------------------------------------

library(tidyverse)
library(shinythemes)
library(shiny)

shiny_biomart_table <- readRDS("data/shiny_biomart_table.rds")


# Define UI for data upload app -------------------------------------------

ui <- fluidPage(theme = shinytheme("flatly"),

    # App title
    titlePanel(h1("Human Gene ID Converter")),

    # Sidebar layout with input and output definitions
    sidebarLayout(

        # Sidebar panel for inputs
        sidebarPanel(

            # Input: Select a file
            fileInput("file1", "Choose text file:",
                      multiple = TRUE,
                      accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")),

            tags$hr(),

           # Input: Select separator ----
            radioButtons("sep", "Separator",
                         choices = c(Comma = ",",
                                     Semicolon = ";",
                                     Tab = "\t"),
                         selected = ","),


           # Input: Checkbox if file has header
           checkboxInput("header", "File contains header", FALSE),

           tags$hr(),

           downloadButton("Matched_dl", "Download Matched Genes"),

           tags$br(),
           tags$br(),

           downloadButton("NoMatch_dl", "Download Non-Matching Genes"),

           tags$br(),
           tags$br(),

           downloadButton("LOC_dl", "Download LOC Genes")

        ),

        # Main panel for displaying outputs
        mainPanel(

            # Output: Data file
            h3("Matching Genes:\n"),
            tableOutput("Matched"),

            h3("Non-matching Genes:\n"),
            tableOutput("NoMatch"),

            h3("LOC Genes:\n"),
            tableOutput("LOC")

        )

    )
)


# Define server logic to read selected file -------------------------------

server <- function(input, output) {


    # Matching Genes ------------------------------------------------------

    Matched <- reactive({

        req(input$file1)

        input_genes <- read.csv(input$file1$datapath,
                                header = input$header,
                                sep = input$sep)


        # Remove any spaces, remove IDs containing a "."
        clean_genes_1 <- as.character(input_genes[,1]) %>%
            str_replace_all(., pattern = "\\s", replacement = "") %>%
            str_subset(., pattern = "\\.", negate = T)

        clean_genes_2 <- grep(clean_genes_1, pattern = "^LOC[0-9]+$", value = T, invert = T)


        # Find matching entries
        matching_genes <- filter_all(shiny_biomart_table, any_vars(. %in% clean_genes_2)) %>%
            distinct(hgnc_symbol, .keep_all = T)


        # Create and return output
        return(matching_genes)
    })


    output$Matched <- renderTable({
        Matched()
    })

    output$Matched_dl <- downloadHandler(
        filename = "matched_genes.csv",
        content = function(file) {
            write.csv(Matched(), file, row.names = FALSE, quote = FALSE)
        }
    )


    # Non-matching Genes --------------------------------------------------

    NoMatch <- reactive({
        req(input$file1)

        input_genes <- read.csv(input$file1$datapath,
                                header = input$header,
                                sep = input$sep)


        # Remove any spaces, remove IDs containing a "."
        clean_genes_1 <- as.character(input_genes[,1]) %>%
            str_replace_all(., pattern = "\\s", replacement = "") %>%
            str_subset(., pattern = "\\.", negate = T)

        clean_genes_2 <- grep(clean_genes_1, pattern = "^LOC[0-9]+$", value = T, invert = T)


        # Find matching entries
        matching_genes <- filter_all(shiny_biomart_table, any_vars(. %in% clean_genes_2)) %>%
            distinct(hgnc_symbol, .keep_all = T)

        # Find input genes with no matches
        matched_chr <- unlist(matching_genes) %>% as.character()
        nonmatch_genes <- data.frame(
            Genes = setdiff(clean_genes_2, matched_chr)
        )

        # Create and return output
        return(nonmatch_genes)
    })


    output$NoMatch <- renderTable({
        NoMatch()
    })

    output$NoMatch_dl <- downloadHandler(
        filename = "non-matching_genes.csv",
        content = function(file) {
            write.csv(NoMatch(), file, row.names = FALSE, quote = FALSE)
        }
    )


    # LOC genes -----------------------------------------------------------

    LOC <- reactive({
        req(input$file1)

        input_genes <- read.csv(input$file1$datapath,
                                header = input$header,
                                sep = input$sep)


        # Remove any spaces, remove IDs containing a "."
        clean_genes_1 <- as.character(input_genes[,1]) %>%
            str_replace_all(., pattern = "\\s", replacement = "") %>%
            str_subset(., pattern = "\\.", negate = T)


        # Get the LOC IDs
        LOC_genes <- data.frame(
            Genes = grep(clean_genes_1, pattern = "^LOC[0-9]+$", value = T)
        )

        # Create and return output
        return(LOC_genes)
    })

    output$LOC <- renderTable({
        LOC()
    })

    output$LOC_dl <- downloadHandler(
        filename = "LOC_genes.csv",
        content = function(file) {
            write.csv(LOC(), file, row.names = FALSE, quote = FALSE)
        }
    )

}


# Run the app -------------------------------------------------------------

shinyApp(ui, server)
