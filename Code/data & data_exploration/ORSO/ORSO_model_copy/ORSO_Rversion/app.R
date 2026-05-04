library(shiny)

runApp(list(
  ui=pageWithSidebar(headerPanel("Introductions"),
                     sidebarPanel(textInput("text1", "Column 1"),
                                  numericInput("text2", "Column 2", 1, min = 1, max = 100),
                                  actionButton("update", "Update Table"),
                                  actionButton("clear", "Clear Table")),
                     mainPanel(tableOutput("table1"))),
  server=function(input, output, session) {
    values <- reactiveValues()
    values$df <- data.frame(Column1 = numeric(0), Column2 = numeric(0))
    ModTable1 <- observe({
      if(input$update > 0) {
        newLine <- isolate(c(input$text1, input$text2))
        isolate(values$df[nrow(values$df) + 1,] <- c(input$text1, input$text2))
      }
    })
    ModTable2 <- observe({
      if(input$clear > 0) {
        isolate(values$df <- data.frame(Column1 = numeric(0), Column2 = numeric(0)))
      }
    })
    output$table1 <- renderTable({values$df})
  }))