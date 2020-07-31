

library(shiny)

# Define UI for application that draws a histogram
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GEOquery")

library(GEOquery) 
library(tidyverse)
library(stringr)

all_data = getGEO("GSE21935") 
all_data = all_data$GSE21935_series_matrix.txt.gz 

pheno = all_data@phenoData@data 
assay = assayDataElement(all_data@assayData, 'exprs') 


x <- as.data.frame(t(head(assay, 10))) 


# extract names of genes
gene_names <- colnames(x)

ui <- fluidPage(
    # Application title
    titlePanel("Visualisation of Gene expression data"),
    
    # Sidebar layout with input and output definitions ---
    sidebarLayout(
        # Sidebar with a slider and selection inputs
        sidebarPanel(
            selectInput("selection", "Choose a gene:",
                        choices = gene_names),
        ),
        
        # Main panel for displaying output
        mainPanel(
            # Output: Density plots
            plotOutput(outputId = "plot")
        )
    )
)




y <-  rep("NA", dim(x)[1])

# creating a for loop to assign 0 to control and 1 to scz
for (i in 1:(dim(x)[1])) {
    if(pheno$characteristics_ch1.4[i]=="disease state: control"){
        y[i] = 0
    } else {
        y[i] = 1
    }  
}

# combining data from x and y so that it is easier to specify conditions
z <- cbind(x,y)

# Define server logic required to generate plots ----
server <- function(input, output) {
    
    kde_plots <- function(selection){
        gene_num <- grep(selection, gene_names)
        kde_left <- density(z[z$y == "0",gene_num])
        kde_right <- density(z[z$y == "1",gene_num])
        return(list(kde_left,kde_right))
    }
    
    output$plot <- renderPlot({
        par(mfrow=c(2,1))
        plot(kde_plots(input$selection)[[1]], main = "Kernel density plot for control group")
        plot(kde_plots(input$selection)[[2]], main = "Kernel density plot for schizophrenic group")
    })
}


# Run the application 
shinyApp(ui = ui, server = server)

