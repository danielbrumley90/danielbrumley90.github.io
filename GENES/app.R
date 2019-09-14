############################## PACKAGES ##############################

library(shiny)
library(randomForest)
library(shinyRGL)
library(rglwidget)
library(forestFloor)
library(dplyr)
library(purrr)
library(readr)
library(stringr)

############################## DATA/RESULTS ##############################

# pathway data
pathways <- read_file("http://people.duke.edu/~hp44/data/breast.gmt") %>% # read data into string
  str_replace_all("\n", "") %>% # remove new line characters
  str_replace_all("\tnone", "") %>% # remove "none" genes
  str_split("\t\r") %>% # split string into a list of pathways
  map(str_replace_all, "\t", ", ") %>% # replace all tabs
  map(str_split_fixed, ", ", 2) %>%  # split pathway names and genes
  data.frame(stringsAsFactors = FALSE) %>% # store as a data frame
  head(-1) # remove the last row (empty)

colnames(pathways) <- c("Pathway", "Genes")

# gene expression data
gene.expression <- read_delim("http://people.duke.edu/~hp44/data/breast.gct", "\t", skip = 2)
gene.names <- gene.expression$probeID
gene.expression <- data.frame(t(gene.expression[, -c(1, 2)]))
colnames(gene.expression) <- gene.names

# response
response <- read_file("http://people.duke.edu/~hp44/breast.cls") %>%
  str_extract("(\\d\\s){48}\\d$") %>%
  str_split(" ", simplify = TRUE) %>%
  factor(labels = c("AT", "BT", "LT"))

# store glycolysis and hem pathway information
glycolysis <- dplyr::select(gene.expression, one_of(str_split(pathways$Genes[pathways$Pathway == "Glycolysis - Glucone"], ", ") %>% unlist))
hem <- dplyr::select(gene.expression, one_of(str_split(pathways$Genes[pathways$Pathway == "BC-Regulation of hem"], ", ") %>% unlist))

# suppress errors
options(shiny.sanitize.errors = TRUE)

# store pathways
pathways <- list(glycolysis, hem)

# store response
y <- response

# store random forest results
rf.results <- list(randomForest(pathways[[1]], y, keep.forest = TRUE, keep.inbag = TRUE, samp = 20, importance = TRUE),
                   randomForest(pathways[[2]], y, keep.forest = TRUE, keep.inbag = TRUE, samp = 20, importance = TRUE))

# store forest floor results
ff.results <- list(forestFloor(rf.results[[1]], pathways[[1]]),
                   forestFloor(rf.results[[2]], pathways[[2]]))

############################## USER INTERFACE ##############################

ui <- fluidPage(
  
  headerPanel("Random Forest Visualization with Breast Cancer Data"),
  
  sidebarPanel(
    selectInput("pathway", 
                label = "Choose Pathway:", 
                choices = c("Glycolysis-Gluconeogenesis" = "1", "BC-Regulation of hem" = "2")),
    uiOutput("geneControls")
  ),
  
  mainPanel(
    tags$style(type="text/css",".shiny-output-error { visibility: hidden; }",".shiny-output-error:before { visibility: hidden; }"),
    tabsetPanel("tabs",
      tabPanel("3D Visualization", value = "3d", webGLOutput("plot3d", width = 800, height = 800)),
      tabPanel("2D Visualization", value = "2d", plotOutput("plot2d"))
    )
  )
)

############################## SERVER LOGIC ##############################

server <- function(input, output, session) {
  
  # selected pathway
  pathway <- reactive({
    pathways[[as.numeric(input$pathway)]]
  })
  
  # random forest results (depends on pathway)
  rf <- reactive({
    rf.results[[as.numeric(input$pathway)]]
  })
  
  # forestFloor results (depends on pathway)
  ff <- reactive({
    ff.results[[as.numeric(input$pathway)]]
  })
  
  output$geneControls <- renderUI({
    if (input$tabs == "3d") {
      return(selectInput("gene1", 
                         label = "Select Gene 1:", 
                         choices = NULL))
    } else {
      return(selectInput("gene2", 
                         label = "Select Gene 2:", 
                         choices = NULL))
    }
  })

  # update gene1 and gene2 selectInput choices (depends on pathway)
  observe({
    
    # get current pathway
    pathway <- pathway()
    
    # update gene1 selectInput
    updateSelectInput(session, "gene1", 
                      choices = names(pathway), 
                      selected = names(pathway)[1])
    
    # update gene2 selectInput
    updateSelectInput(session, "gene2", 
                      choices = names(pathway), 
                      selected = names(pathway)[2])
  })
  
  # selected genes (depends on gene1, gene2)
  genes <- reactive({
    pathway()[, c(input$gene1, input$gene2)]
  })
  
  # indices of selected genes (depends on gene1, gene2)
  gene.indices <- reactive({
    c(which(names(pathway()) == input$gene1),
      which(names(pathway()) == input$gene2))
  })
  
  # render 3D plot
  output$plot3d <- renderWebGL({
    show3d(ff(), gene.indices())
  })
  
  # render 2D plot
  output$plot2d <- renderPlot({
    plot(genes(), col = rf()$predicted, pch = 20, cex = 3)
  })
}

############################## RUN APPLICATION ##############################

shinyApp(ui = ui, server = server)