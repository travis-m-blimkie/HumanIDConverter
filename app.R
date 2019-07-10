
# Load libraries, global variables ----------------------------------------

library(tidyverse)
library(shiny)

shiny_biomart_table <- readRDS("data/shiny_biomart_table.rds")


# Define UI for data upload app -------------------------------------------
ui <- fluidPage(

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

           tags$hr(),

           # Input: Checkbox if file has header
           checkboxInput("header", "File contains header", FALSE)


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


    output$Matched <- renderTable({


        req(input$file1)

        input_genes <- read.csv(input$file1$datapath,
                                header = input$header,
                                sep = input$sep)


        # Remove any spaces, remove IDs containing a "."
        clean_genes_1 <- as.character(input_genes[,1]) %>%
            str_replace_all(., pattern = "\\s", replacement = "") %>%
            str_subset(., pattern = "\\.", negate = T)

        clean_genes_2 <- grep(clean_genes_1, pattern = "^LOC", value = T, invert = T)


        # Find matching entries
        matching_genes <- filter_all(shiny_biomart_table, any_vars(. %in% clean_genes_2)) %>%
            distinct(hgnc_symbol, .keep_all = T)


        # Create and return output
        return(matching_genes)

    })

    output$NoMatch <- renderTable({


        req(input$file1)

        input_genes <- read.csv(input$file1$datapath,
                                header = input$header,
                                sep = input$sep)


        # Remove any spaces, remove IDs containing a "."
        clean_genes_1 <- as.character(input_genes[,1]) %>%
            str_replace_all(., pattern = "\\s", replacement = "") %>%
            str_subset(., pattern = "\\.", negate = T)

        clean_genes_2 <- grep(clean_genes_1, pattern = "^LOC", value = T, invert = T)


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


    output$LOC <- renderTable({


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
            Genes = grep(clean_genes_1, pattern = "^LOC", value = T)
        )

        # Create and return output
        return(LOC_genes)

    })

}


# Run the app -------------------------------------------------------------
shinyApp(ui, server)
