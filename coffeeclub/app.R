library(shiny)
library(RMySQL)
library(DBI)

shinyApp(
  ui = shinyUI(
      fluidPage(theme ="darkly.css",
      titlePanel("Coffee Club Data"),
      sidebarLayout(
        sidebarPanel(
                fluidRow(
                      h4(uiOutput("userID")),
                      h4(dateInput("date", "Date:", value = NULL, min = NULL, max = NULL, startview = "month", weekstart = 0,language = "en", width = NULL)),
                      h4(radioButtons("coffee", "Coffee:", c( "Single" = 1,"Double" = 2), selected="1", inline = TRUE)),
                      br(),
                      #h3(actionButton("enterdata", "Enter", style='font-size:100%')),
                      actionButton("enterdata", "", icon = icon("play-circle")),
                      tags$button(id = "singlecoffee", class = "btn action_button",tags$img(src = "coffee.png",height = "50px"))
                        )
                      ),
        mainPanel(
                fluidRow(
                      br(),
                      HTML("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"),
                      h5(textOutput("text")),
                      h5(textOutput("text1")), h5(tableOutput("tbl")), h5(tableOutput("tbl1"))
                        )
                  )
                )
              )
          ),

  
  server = function(input, output, session){
    out <- reactiveVal("Tap Coffee Button !!!!")
    
    
    
    data_sets  <- reactive({
      conn1 <- dbConnect(
        drv = RMySQL::MySQL(),
        dbname = "coffeeclub",
        host = "localhost",
        username = "root",
        password = "abc123")
      on.exit(dbDisconnect(conn1), add = TRUE)
      
      return (dbGetQuery(conn1, paste0("SELECT  DISTINCT DisplayName  FROM coffeeclub.coffeeclub_members  ORDER BY DisplayName ASC;")))
      
    })
    
    # Drop-down selection box for which data set
    output$userID <- renderUI({
      selectInput("uID", "User ID", as.list(data_sets()))
    })
    
    output$data_table <- renderTable({
      get(input$uID)
    })
    
    
    goButton1 <- eventReactive(input$enterdata,{
      output$tbl <- renderTable({
        conn <- dbConnect(
          drv = RMySQL::MySQL(),
          dbname = "coffeeclub",
          host = "localhost",
          username = "root",
          password = "abc123")
        on.exit(dbDisconnect(conn), add = TRUE)
        dbGetQuery(conn, paste0("SELECT UserID as Name, date as Date, coffee as Coffee FROM coffeeclub WHERE date >= CURDATE()  ORDER BY ID DESC LIMIT 10;"))
      })
    })
    
    goButton2 <- eventReactive(input$enterdata,{
      
      conn1 <- dbConnect(
        drv = RMySQL::MySQL(),
        dbname = "coffeeclub",
        host = "localhost",
        username = "root",
        password = "abc123")
      on.exit(dbDisconnect(conn1), add = TRUE)
      dbGetQuery(conn1, paste0("insert into coffeeclub ( UserID, date, coffee) values ('",input$uID,"', '",input$date,"',",input$coffee, ") ;"))
      
    })
    
    output$text1 <- renderText({
      goButton1()
      goButton2()
    })
    
  
  }
)