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
ui <- fluidPage(

   # Application title
   titlePanel("Analýza citlivosti povodí na klimatickou změnu"),

   # Sidebar with a slider input for number of bins
   sidebarLayout(
      sidebarPanel(

      ),

      # Show a plot of the generated distribution
      mainPanel(
        plotOutput(
          'mean-sd',
          brush = brushOpts(
            id = 'Brush1'
          ),
          dblclick = dblclickOpts(id = 'DblClick1'),
          click = 'Click1'
        )
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  output[['mean-sd']] <- renderPlot({

    ggplot(vstat) + geom_point(aes(x = P, y = T, shape = EXP, col = PER), size = 3, alpha = .5) + facet_wrap(~variable, scale = 'free', ncol = 2) + if (is.null(input$Click1)) {NULL} else {geom_point(aes(x = input$Click1$x, y = input$Click1$y), size = 4)}
  })
}

# Run the application
shinyApp(ui = ui, server = server)

