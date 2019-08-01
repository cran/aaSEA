ui <- dashboardPage(
  dashboardHeader(title = "Substitution Effect Analysis"),
# begin of dashboard sidebar contents-----
    dashboardSidebar(
    # sidebar menu contents -----
    sidebarMenu(
                menuItem(text = "START", tabName = "begin", icon = icon('check')),
                menuItem(text = 'Substitution, Property Change', tabName = 'pc', icon = icon('cog')),
                menuItem(text = 'Coevolving Site Computation', tabName = 'cm', icon = icon('cog')),
                menuItem(text = 'Coevolving Property Change', tabName = 'cmpc', icon = icon('cog'))
               )
),
# begin of dashboard body desing-----
  dashboardBody(
    tabItems(
     # first tab content-----
       tabItem(tabName = "begin",
              # output About text box
               uiOutput("abtxt"),
              # output file uploader
              box(title = tags$b("Upload alignment"), solidHeader = TRUE, status = "info",
                  fileInput(inputId = "ali",
                            label =  "Choose Alignment file",
                            multiple = F, accept = ".fasta"))
              ),
      # Second tab content regarding property changes-----
      tabItem(tabName = 'pc',
              fluidRow(
              box(title = "Choose sites with", solidHeader = TRUE,status = "success",
                  width = 3, height = '150px',
                  radioButtons("Sel", NULL,
                               choices = list("Single Substitution" = "ss",
                                              "Multiple Substitutions" = "ms"), selected = "ss")
              ),
              # output choices box for peoperty changes
              box(title = "Get changes in", solidHeader = TRUE, status = 'primary', width = 3, hight = "150px", collapsed = FALSE, collapsible = TRUE,
              radioButtons("SelChange", NULL,
                           choices = list("Cruciani","Kidera", "Fasgai", "AAindex"
                                          ), selected = "Cruciani")
              ),
              
              box(title = "Select Property", solidHeader = TRUE, status = 'success', width = 6, hight = "100px", collapsed = FALSE, collapsible = TRUE, 
                  # selectInput("PropIndex", label = "Rownames", choices = "rows")
                  uiOutput('rows')
                  )),
              fluidRow(
              box(title = "Property Changes", width = 12,height = '500px', collapsible = TRUE, collapsed = FALSE, solidHeader = TRUE, status = 'success',
              withSpinner(plotlyOutput("pcplot"))),
              box(width = 12, withSpinner(uiOutput("table2")))
              )),
    
# Third tab content. Regarding the Correlated mutations------
      tabItem(tabName = 'cm',
              # output the choices box for correlated mutations
              fluidRow(
              box(title = "Correlated mutations based on", solidHeader = TRUE, status = 'primary', width = 4,
                  radioButtons('cSel', NULL,
                              choices = list('McLachlan Based Substitution Correlation' = 'mcbasc',
                                             'Mutual Information Product' = 'mip',
                                             'Observed Minus Expected Squared' = "omes",
                                             'Explicit Likelihood of Subset Covariation' = "elsc"),
                              selected = "mcbasc")
                  )),
              fluidRow(
                box(title = "Simple network diagram",solidHeader = TRUE, width = 6, collapsible = TRUE, collapsed = FALSE, status = "success",height = '500px',
                    withSpinner(simpleNetworkOutput("simple"))
                    ),
                box(withSpinner(uiOutput("table3")))
              )
      ),
# Fouth tab content. Regarding Property changes in correlated mutations-------
     tabItem(tabName = 'cmpc',
             # output choice boxes
             fluidRow(
             box(title = "Get changes in", solidHeader = TRUE, status = 'primary', width = 3, hight = "500px", collapsed = FALSE, collapsible = TRUE,
             radioButtons("CorSelChange", NULL,
                              choices = list("Cruciani","Kidera", "Fasgai", "AAindex"
                              ), selected = "Cruciani")
             ),
             
             box(title = "Select Property", solidHeader = TRUE, status = 'primary', width = 6, hight = "500px", collapsed = FALSE, collapsible = TRUE, 
                 uiOutput('prows')
             )),
               fluidRow(
               box(title = "Correlated Property Changes", width = 12,height = '500px', collapsible = TRUE, collapsed = FALSE, solidHeader = TRUE, status = 'success',
                   withSpinner(plotlyOutput("ccorrPlot")))
             ),
             fluidRow(
               withSpinner(uiOutput("table4")))
    )
    )
    )
    )