#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(theme ="bootstrap.css",
                titlePanel("Please enter data in this web-form !! One Entry per project"),
                sidebarLayout(
                  sidebarPanel(
                    selectInput("group", "Choose a Group:", c("UG"= "UG", "CEL" = "CEL"), selected=NULL),
                    selectInput("member", "Choose a user:", c("Vimal"= "Vimal", "Alex" = "Alex"), selected=NULL),
                    actionButton("showU","Show Data"), 
                    br(),br(),
                   
                    h4("Last 10 entries!!"), 
                    textOutput("text10"),
                    tableOutput("viewLast")
                    ,width =2),
                  
                  mainPanel(tabsetPanel( 
                                         tabPanel(h3("New Users"),
                                                  h3(selectInput("project", "Choose Project:",c("none" = "none","URPP" = "URPP","ERC" = "ERC"), selected=NULL)),
                                                  
                                                  h4(textInput("shortname", "Title (eg. MEDIA, TRANSNET)","")),
                                                  textAreaInput(inputId="fullname", "Full Name (eg. Transcriptional Network analysis of Developing Female Gametophyte)",rows=3, width='300%', resize = "both", value=""),
                                                  
                                                  h4(textAreaInput(inputId="abstract", "Abstract:",rows=20, width='300%', resize = "both", value="")),
                                                  
                                                  h4(textInput("datatype", "Sequencing data Type (eg. RNA_seq)","")),
                                                  h4(textInput("dept", "Department (eg. UG)","")),
                                                  h4(textInput("phone", "Phone Number","+41 ")),
                                                  h4(numericInput("age", "Age", 30, min = 18, max = 110)),
                                                  h4(textInput("jobtitle", "Job title (eg. Postdoc)","")),
                                                  h4(textInput("nationality", "Nationality (eg. German)","")),
                                                  h4(dateInput("membersince", "MemberSince (eg. DATE)")),
                                                  
                                                  actionButton("addU","Add member")
                                                  , width=4)
                                         
                                         
                  ), width=6)
                )
)





# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      x    <- faithful[, 2] 
      bins <- seq(min(x), max(x), length.out = input$bins + 1)
      
      # draw the histogram with the specified number of bins
      hist(x, breaks = bins, col = 'darkgray', border = 'white')
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

