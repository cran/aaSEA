server <- function(input, output) {
  storeWarn<- getOption("warn")
  options(warn = -1) 
  
  #about text output
  output$abtxt <- renderUI({
    box(
      title = tags$b("About aaSEA"),
      solidHeader = TRUE,
      status = "warning",
      "* Upload FASTA formated Multiple Sequence Alignment ",
      br(),
      "* Use sidebar menu to navigate specific analysis",
      br(),
      "* Use filters on respective tabs to change outputs"
    )
  })
  #################################
  ### Prepare substitution Data ###
  #################################
  # Single vs multiple substitution capture
  ssda <- reactive({
    req(input$ali$datapath, cancelOutput = FALSE)
    mut <- getAASub(input$ali$datapath)
    if (input$Sel == "ss") {
      sf <- mut$singleSub
    } else if (input$Sel == "ms") {
      sf <- mut$multiSub
    } 
    return(sf)
  })
  #  Physicochemical property changes----
  ## Filters for property selection----
  output$rows <- renderUI({
    mydata <- readRDS(paste0("data/", input$SelChange, ".rds"))
    selectInput(inputId = "PropIndex",
                label = "Property in Property Group",
                choices = rownames(mydata))
  })
  #Physico chemical property change claculation
  chaData <- reactive({
    myD <- ssda()
    propDF <- readRDS(paste0("data/", input$SelChange, ".rds", sep = ""))
    pDFname <- input$SelChange
    myIndex <- input$PropIndex
    pIndex <- match(myIndex, rownames(propDF))
    changes <- 
      getPropChange(
        subFile = myD,
        propertyDF = pDFname,
        propertyIndex = pIndex
      )
    return(changes)
  })
  ################################
  ### Correlated mutation data ###
  ################################
  # calculate site correlations
  corsubDat <- reactive({
    req(input$ali$datapath, cancelOutput = FALSE)
    myCorMe <- input$cSel # get method to calculate correlations
    selMat <- getCorSites(fileLoc = input$ali$datapath, corMethod = myCorMe) # select columns based on method selected
    res <- getTopSub(selMat = selMat)
    return(list("SelMatrix" = selMat,
                "topSub" = res))
  })
  # Property changes in correlated mutations----
  # Filters for correlated sites
  output$prows <- renderUI({
    mydata <- readRDS(paste0("data/", input$CorSelChange, ".rds"))
    selectInput(inputId = "corPropIndex",
                label = "Property in Property Group",
                choices = rownames(mydata))
  })
  corCha <- reactive({
    mycorDa <- corsubDat()$topSub
    head(mycorDa)
    cpropDF <- readRDS(paste0("data/", input$CorSelChange, ".rds", sep = "")) # read data for row names
    cpDFname <- input$CorSelChange # Just name of the data frame # not working
    cmyIndex <- input$corPropIndex # not working
    cpIndex <- match(cmyIndex, rownames(cpropDF)) # not working
    corChanges <- getCorPropChange(corSubFile = mycorDa, propertyDF = cpDFname,propertyIndex = cpIndex)
    return(corChanges)
  })
  # Try to display only the property correlations with selected sites
  output$prows <- renderUI({
    mydata <- readRDS(paste0("data/", input$CorSelChange, ".rds"))
    selectInput(inputId = "corPropIndex",
                label = "Property in Property Group",
                choices = rownames(mydata))
  })
  propCor <- reactive({
    req(input$ali$datapath, cancelOutput = FALSE)
    myCorMe <- input$cSel
    myselMat <- corsubDat()$SelMatrix
    #selMat <- getCorSites(fileLoc = input$ali$datapath, corMethod = myCorMe)
    cpropDF <- readRDS(paste0("data/", input$CorSelChange, ".rds", sep = "")) # read data for row names
    cpDFname <- input$CorSelChange # Just name of the data frame # not working
    cmyIndex <- input$corPropIndex # not working
    cpIndex <- match(cmyIndex, rownames(cpropDF)) # not working
    corDF <-
      getPropCorr(selMat = myselMat,
                  propertyDF = cpDFname,
                  propertyIndex = cpIndex)
    return(corDF)
  })
  ##########################
  ### output definitions ###
  ##########################
  # property changes-----
  output$table2 <- renderUI({
   # shiny::validate(need(nrow(chaData()) != 0, message = "No data to display. upload alignment"))
    output$propCha <- DT::renderDataTable(
      chaData(),
      filter = "top",
      rownames = FALSE,
      options = list(
        lengthMenu = c(5, 10, 30, 50),
        pageLength = 5
      )
    )
    box(
      title = paste("Property Changes in ", input$SelChange, " - ", input$PropIndex),
      width = 12,
      collapsible = TRUE,
      collapsed = FALSE,
      status = 'success',
      solidHeader = TRUE,
      div(style = 'overflow-x:scroll',
          fluidRow(
            column(
              width = 12,
              downloadButton(outputId = "propChange", label = "Download"),
              DT::dataTableOutput("propCha")
            )
          ))
    )
  })
  output$propChange <-
    downloadHandler(
      filename = function() {
        paste(input$Sel, "_", input$SelChange, ".csv", sep = "")
      },
      content = function(file) {
        write.csv(chaData(), file, row.names = FALSE)
      }
    )
  # Plot for property changes
  output$pcplot <- renderPlotly({
    mut <- getAASub(input$ali$datapath)
    # if single substitution is selected
    if (input$Sel == "ss") {
      sdf <- mut$singleSub
        propDF <-
          readRDS(paste0("data/", input$SelChange, ".rds", sep = ""))
        pDFname <- input$SelChange
        myIndex <- input$PropIndex
        pIndex <- match(myIndex, rownames(propDF))
        sdfchanges <-
          getPropChange(
            subFile = sdf,
            propertyDF = pDFname,
            propertyIndex = pIndex
          )
        
        ss <- sdfchanges
        ss <- subset(ss, ss$wt != "-")
        ss <- subset(ss, ss$mu != "-")
        sso <- ss[order(ss$Delta.Prop), ]
        sso$substitution <- factor(sso$substitution,
                                   levels = unique(sso$substitution)[order(sso$Delta.Prop,
                                                                           decreasing = TRUE)])
        m <-
          list(l = 50,
               r = 20,
               b = 50,
               t = 20) # l = left; r = right; t = top; b = bottom
        #mycolors <- colorRampPalette(RColorBrewer::brewer.pal(8,"Set2"))(~mu)
        plot_ly(sso,x = ~ substitution,y = ~ Delta.Prop,type = "bar",color = ~ mu) %>%
        #plot_ly(sso,x = ~ substitution,y = ~ Delta.Prop,type = "bar",color = mycolors) %>%
          layout(xaxis = list(tickangle = 45),
                 margin = m,
                 title = "Single substitutions and associated property changes")
    } else if (input$Sel == "ms") {
      # if Multiple seubstitutions are selected
      ms <- mut$multiSub
      propDF <-
        readRDS(paste0("data/", input$SelChange, ".rds", sep = ""))
      pDFname <- input$SelChange
      myIndex <- input$PropIndex
      pIndex <- match(myIndex, rownames(propDF))
      print(ms)
      msc <-
        getPropChange(
          subFile = ms,
          propertyDF = pDFname,
          propertyIndex = pIndex
        )
      print(msc)
      msco <- subset(msc, msc$wt != "-")
      msco <- subset(msco, msco$mu != "-")
      msco$xlab <- paste(msco$wt, msco$site, sep = "")
      
      # simple Heat map route
      plot_ly(
        x = as.character(msco$xlab),
        y = msco$mu,
        z = round(msc$Delta.Prop, 2),
        type = "heatmap",
        zauto = FALSE
      ) %>%
        layout(title = "Multiple substitutions and associated property changes")
    }
  })
  #================================
  # correlated mutations-----
  #================================
  output$table3 <- renderUI({
   # shiny::validate(need(nrow(corsubDat()$topSub) != 0, message = "No data to display. upload alignment"))
    output$corSub <- DT::renderDataTable(
      corsubDat()$topSub,
      filter = "top",
      rownames = FALSE,
      options = list(
        lengthMenu = c(5, 10, 30, 50),
        pageLength = 5
      )
    )
    box(
      title = "Co-evolving site pairs",
      width = 12,
      height = '500px',
      collapsible = TRUE,
      collapsed = FALSE,
      solidHeader = TRUE,
      status = 'success',
      div(style = 'overflow-x:scroll',
          fluidRow(
            column(
              width = 12,
              downloadButton(outputId = "corsubsites", label = "Download"),
              DT::dataTableOutput("corSub")
            )
          ))
    )
  })

  
  
  output$corsubsites <-
    downloadHandler(
      filename = function() {
        paste(input$cSel, ".csv", sep = "")
      },
      content = function(file) {
        write.csv(corsubDat()$topSub, file, row.names = FALSE)
      }
    )
  output$simple <- renderSimpleNetwork({
    plotDF <- corsubDat()$topSub
    simpleNetwork(plotDF, fontSize = 12, fontFamily = "sans-serif", zoom = TRUE)
  })
  
  # correlated mutation property changes-----
  output$table4 <- renderUI({
    #shiny::validate(need(nrow(corCha()) != 0, message = "No data to display. upload alignment"))
    output$corPropCha <- DT::renderDataTable(
      corCha(),
      filter = "top",
      rownames = FALSE,
      options = list(lengthMenu = c(10, 30, 50),
                     pageLength = 10)
    )
    box(
      title = "Correlated substitution Property Changes",
      width = 12,
      collapsible = TRUE,
      collapsed = FALSE,
      solidHeader = TRUE,
      status = 'success',
      div(style = 'overflow-x:scroll',
          fluidRow(
            column(
              width = 12,
              downloadButton(outputId = "corsubproperty", label = "Download"),
              DT::dataTableOutput("corPropCha")
            )
          ))
    )
  })
  output$corsubproperty <-
    downloadHandler(
      filename = function() {
        paste(input$selPc, ".csv", sep = "")
      },
      content = function(file) {
        write.csv(corCha(), file, row.names = FALSE)
      }
    )
  output$ccorrPlot <- renderPlotly({
    mydf <- propCor()
    plot_ly(
      x = as.character(mydf$Pos1),
      y = as.character(mydf$Pos2),
      z = round(mydf$cor, 2),
      type = "heatmap",
      zauto = FALSE
    ) %>%
      layout(title = "Site correlations with respect to Selected property")
  })
}