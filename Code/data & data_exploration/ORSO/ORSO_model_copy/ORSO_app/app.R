#
# This is a Shiny web application for running sea otter population simulations
#  along coastal Oregon
library(shiny)
library(shinythemes)
library(DT)
library(shinyWidgets)
library(shinyBS)
library(shinyjs)
library(rhandsontable)
library(gtools)
library(stats)
library(mvtnorm)
library(boot)
library(ggplot2)
library(mapproj)
library(rgdal)
library(dplyr) 
library(ggrepel)
library(ggsn)
library(reshape2)
#
Bdat = read.csv("./data/OR_pops.csv", header = TRUE)  # Data for coastal blocks
NBlks <- nrow(Bdat)

instructfile = "./data/Instruct.txt"
Blurb1 = readChar(instructfile, file.info(instructfile)$size)
blurbfile = "./data/Blurb.txt"
Blurb2 = readChar(blurbfile, file.info(blurbfile)$size)

Csections = Bdat$Segment[order(Bdat$Ycoord)]
arealist <- as.list(Csections); names(arealist) <- Csections
# User interface part of app
ui <- shinyUI(
    navbarPage(
        theme = shinythemes::shinytheme("sandstone"),
        "Oregon Sea Otter Population Model (ORSO V 1.0)",
        tabPanel("Setup Model Simulations",
                 # Sidebar panel: user input
                 useShinyjs(),
                 sidebarPanel(
                     helpText(Blurb1),
                     div(style="display:inline-block; horizontal-align: left",
                         actionButton("RunSim", "Run Simulations Now", 
                         class = "btn-primary"),width=10),
                     div(style="display:inline-block; margin-left:100px"),
                     div(style="display:inline-block; horizontal-align: right",
                         downloadButton('GetManual', 'Download User Manual')),
                     div(style="margin-top:10px") ,
                     selectizeInput(
                         "InitSect",
                         label = ("Select a coastal section for re-introducing otters (refer to map):"),
                         choices = arealist, multiple = FALSE, width='100%', selected = NULL,# 
                         options = list(placeholder = 'click to select one or more sections',
                                        onInitialize = I('function() { this.setValue(""); }'),            
                         'plugins' = list('remove_button'),'create' = TRUE,'persist' = FALSE)),
                     sliderInput("InitN", "Number of otters initially added to coastal section",
                                 value = 1, min = 1, max = 100, step = 1),
                     bsPopover("InitN",
                               "Adjust numbers of otters: more information",
                               paste0("After selecting a coastal section from a drop-down list (see map for locatations ",
                                      "of each section), use this slider to adjust numbers of otters to be added to ",
                                      "that section in the initial translocation event."),
                               placement = "top", trigger = "hover"),
                     sliderInput("AddN", "Number of otters (per year) in supplemental introductions",
                                 value = 0, min = 0, max = 20, step = 1),
                     bsPopover("AddN",
                               "Adjust supplemental otter introductions: more information",
                               paste0("After selecting a coastal section from a drop-down list (see map for locatations ",
                                      "of each section) and initial translocation numbers, use this slider to adjust ",
                                      "annual number of otters to be added to this section in subsequent years."),
                               placement = "top", trigger = "hover"),                     
                     actionButton("update", "Update Introduction Table"),
                     bsPopover("update",
                               "Update table: more information",
                               paste0("After selecting a coastal section and numbers of otters to be introducted, click this button ",
                                      "to add these translocation parameters to the simulation. If desired, additional sections ",
                                      "can then be added to the table as secondary translocation sites. "),
                               placement = "top", trigger = "hover"),
                     actionButton("clear", "Clear Introduction Table"),
                     bsPopover("clear",
                               "Clear table: more information",
                               paste0("Use this button to clear translocation parameters from the table and start again. "),
                               placement = "top", trigger = "hover"),                     
                     tableOutput("Intro_table"),
                     div(style="margin-top:10px") ,
                     uiOutput("Reps_slider"),
                     bsPopover("Reps_slider",
                               "Number of iterations: more information",
                                paste0("Increasing the number of replications of a simulation improves the precision ",
                                       "of model predictions, but will take longer to run. At least 100 reps is suggested"),
                               placement = "top", trigger = "hover"),
                     uiOutput("Nyrs_slider"), 
                     uiOutput("Addyrs_slider"),
                     uiOutput("PpnFem_slider"),
                     uiOutput("PpnAd_slider"),
                     uiOutput("Ephase_slider"),
                     uiOutput("AddMortE_slider"),
                     uiOutput("ProbMv_slider"),
                     uiOutput("AgeMvAdj_slider"),
                     uiOutput("EstryMvAdj_slider"),
                     uiOutput("Mvmort_slider"),
                     uiOutput("V_sp_slider"),
                     uiOutput("Rmax_slider"),                     
                     uiOutput("Estoch_slider"),
                     uiOutput("Theta_slider"),
                     width = 4
                 ),
                 # Main panel: map of coastal sections 
                 mainPanel(
                     helpText(Blurb2),
                     img(src='OR_map.png', width = 1100,height = 800,
                         align = "center"),
                     width = 8
                 )
        ),
        tabPanel("Model Output GRAPHS",             
                 tabsetPanel(
                         tabPanel("Population Trend",
                                  plotOutput(outputId="sumgraph1",height = 800,width = 1000)),
                         tabPanel("Range Expansion",
                                  plotOutput(outputId="sumgraph2",
                                             height = 800,width = 1000)),
                         tabPanel("Density Map",
                                  plotOutput(outputId="sumgraph3",
                                             height = 800,width = 600))
                     )
                 ),
        tabPanel("Model Output TABLES",             
                 tabsetPanel(
                     tabPanel("Table 1: Projected Sea Otter Abundance by Year",
                              downloadButton('download1',"Download Table 1"),
                              tableOutput('sumtable1')),
                     tabPanel("Table 2: Projected Abundance by Coastal Section in Final Year", 
                              downloadButton('download2',"Download Table 2"),
                              tableOutput('sumtable2'))
                 )
        )
    )
)
# Server part of app
sv <- shinyServer(function(input, output, session){
    # Create reactive values object to store sim results data frame
    values <- reactiveValues(Pop_Overall = NULL,dfDens = NULL,Tab1=NULL,
                             Hab_Blocks = NULL,Hab_Blocks_Fin = NULL)
    values$Intro_df <- data.frame(Intro_Section = numeric(0), Initial_N = numeric(0), 
                                  Supplemental_N = numeric(0))
    source('runsims.r', local = TRUE)
    #
    ModTable1 <- observe({
        if(input$update > 0) {
            newLine <- isolate(c(input$InitSect, input$InitN, input$AddN))
            isolate(values$Intro_df[nrow(values$Intro_df) + 1,] <- c(input$InitSect, input$InitN, 
                                                                     input$AddN))
            updateSelectizeInput(session = getDefaultReactiveDomain(),
                                 inputId = "InitSect",
                                 choices = arealist,  
                                 options = list(placeholder = 'click to select one or more sections',
                                             onInitialize = I('function() { this.setValue(""); }')) )
            updateSliderInput(session = getDefaultReactiveDomain(),
                              "InitN", value=1 )
            updateSliderInput(session = getDefaultReactiveDomain(),
                              "AddN", value=0 )
        }
    })
    ModTable2 <- observe({
        if(input$clear > 0) {
            isolate(values$Intro_df <- data.frame(Intro_Section = numeric(0), Initial_N = numeric(0),
                                                  Supplemental_N = numeric(0)))
        }
    })
    # Allow user manual to be downloaded
    output$GetManual <- downloadHandler(
        filename = "ORSO_Manual.pdf",
        content = function(file) {
            file.copy("www/Manual.pdf", file)
        })
    # Run sims when button pushed (but only enable if coastal sections selected)
    observe({
        toggleState(id = "RunSim", condition = nrow(values$Intro_df)>0)
    })
    observeEvent(input$RunSim, {
        req(nrow(values$Intro_df)>0) # req(input$InitSect)
        datasim <- runsims() 
        Pop_Overall <- values$Pop_Overall 
        dfDens <- values$dfDens        
        Hab_Blocks <- values$Hab_Blocks
        Hab_Blocks_Fin <- values$Hab_Blocks_Fin
        Nyrs <- nrow(Pop_Overall)
        titletxt <- paste0("Projected Sea Otter Population, ", Nyrs," Years")
        subtxt <- paste0("Average expected abundance (line) ",
                            "with 95% CI for the mean (dark shaded band) ",
                            "and 95% CI for projection uncertainty (light shaded band)")
        maxN <- ceiling( max(Pop_Overall$upper)/100)*100
        #maxD <- ceiling(100*max(dfDens$Density))/100
        maxD <- ceiling(100*max(dfDens$Number))/100
        plt1 <- reactive({
            ggplot(Pop_Overall, aes(Year, Mean))+
                geom_line(data=Pop_Overall)+
                geom_ribbon(data=Pop_Overall,aes(ymin=lower,ymax=upper),alpha=0.2)+
                geom_ribbon(data=Pop_Overall,aes(ymin=CImeanLo,ymax=CImeanHi),alpha=0.3)+
                ylim(0,maxN) +  
                xlab("Year") +
                ylab("Expected sea otter abundance") +
                ggtitle(titletxt, subtitle=subtxt) + 
                theme_classic(base_size = 13)
        })
        plt2 <- reactive({
            ggplot(dfDens, aes(Year, Block)) +
                geom_tile(aes(fill = Number), color = "white") +
                scale_fill_gradient(low = "white", high = "steelblue",limits=c(0, maxD)) +
                xlab("Year") +
                ylab("Coastal section") +
                theme(legend.title = element_text(size = 13),
                      legend.text = element_text(size = 13),
                      plot.title = element_text(size=15),
                      axis.title=element_text(size=13),
                      axis.text.y = element_text(size=8),
                      axis.text.x = element_text(angle = 90, hjust = 1)) +
                labs(fill = "Mean Expected Abundance",
                     title= "Projected Sea Otter Numbers by Coastal Section, by Year") 
        })
        output$sumgraph1 <- renderPlot({
            plt1() 
        })
        output$sumgraph2 <- renderPlot({
            plt2() 
        })
        output$sumgraph3 <- renderImage({
            if (is.null(values$Pop_Overall)) {
                return(NULL)
            } else {
                return(list(
                    src = "./www/mapdens.png",
                    contentType = "image/png",
                    width = 600,
                    height = 800,
                    alt = "Map of Oregon showing projected sea otter abundance"
                ))
            }
            }, deleteFile = TRUE)    
    })
    # Output tables
    #
    output$Intro_table <- renderTable({values$Intro_df})
    #
    # Summary table 1 (rangewide sums by year), updated when sims completed
    output$sumtable1 <- renderTable(values$Tab1) 
    output$download1 <- downloadHandler(
        filename = function(){"Table1.csv"}, 
        content = function(fname){
            write.csv(values$Tab1, fname)
        }
    )
    # Summary table 2 (section sums at Final year), updated when sims completed
    output$sumtable2 <- renderTable(values$Hab_Blocks_Fin) 
    output$download2 <- downloadHandler(
        filename = function(){"Table2.csv"}, 
        content = function(fname){
            write.csv(values$Hab_Blocks_Fin, fname)
        }
    )    
    #
    # User input sliders etc. 
    output$Reps_slider <- renderUI({
        sliderInput("Reps","Number of times to iterate this population simulation",
                     min = 50, max = 1000, value = 250, step = 10, round = TRUE)
    })
    #
    output$Nyrs_slider <- renderUI({
        sliderInput("Nyrs","Number of years for each simulation", min = 10, max = 100,
                    value = 25, step = 5, round = TRUE)
    })
    addPopover(session,"Nyrs_slider",
               "Number of years: more information",
               paste0("The population model will be run for 'N' years into the future. Increasing 'N' can provide ", 
                      "insights into conditions farther in the future, but results become less reliable the farther ",
                      "ahead in time the model is projected."),
               placement = "top", trigger = "hover")
    #
    output$Addyrs_slider <- renderUI({
        sliderInput("AddYrs","Number of years of supplementary introductions", min = 0, max = 25,
                    value = 5, step = 1, round = TRUE)
    })
    addPopover(session,"Addyrs_slider",
               "Number years of supplemental introductions: more information",
               paste0("After initial introduction, more otters may be added in subsequent years " , 
                      "in order to improve success of the establishing sub-population(s). These additional otters ",
                      "could be wild otters or juvenile re-habilitated otters from captivity."),
               placement = "top", trigger = "hover")    
    #
    output$PpnFem_slider <- renderUI({
        sliderInput("PpnFem","Proportion of introduced otters that are female", min = .3, max = .8,
                    value = .6, step = .01, round = TRUE)
    })
    addPopover(session,"PpnFem_slider",
               "Proportion female: more information",
               paste0("Sea otter populations often have a slightly higher proportion of females " , 
                      "than males, and for introductions a higher proportion of females can increase ",
                      "the potential for growth, though there must be some adult males for reproduction. "),
               placement = "top", trigger = "hover")    
    #    
    output$PpnAd_slider <- renderUI({
        sliderInput("PpnAd","Proportion of introduced otters that are adult (vs sub-adult)", min = 0, max = 1,
                    value = .25, step = .01, round = TRUE)
    })
    addPopover(session,"PpnAd_slider",
               "Proportion adult: more information",
               paste0("Only adult sea otters produce pups, so introducing adults can hasten reproduction. " , 
                      "However, in past translocations it has been found that sub-adults may be more likely ",
                      "to succesfully 'take' to their new habitat, so more sub-adults may improve success. "),
               placement = "top", trigger = "hover")    
    #    
    output$Ephase_slider <- renderUI({
        sliderInput("Ephase","Number of years expected population establishment phase", 
                    min = 1, max = 20, value = 10, step = 1, round = TRUE)
    })
    addPopover(session,"Ephase_slider",
               "Number years of population establishment: more information",
               paste0("Newly established sea otter populations often experience an initial period ",
                      "of limited growth and range expansion, as the population becomes established. ",
                      "This establishment period has varied from 5-20 years in previous re-introductions."),
               placement = "top", trigger = "hover")
    #
    output$ProbMv_slider <- renderUI({
        sliderInput("ProbMv","Probability of dispersal during establishment phase (adults)", 
                    min = 0, max = 0.95, value = 0.75, step = .01, round = TRUE)
    })
    addPopover(session,"ProbMv_slider",
               "Probability of establishment phase dispersal: more information",
               paste0("In several previous sea otter translocations, a substantial proportion of the ",
                      "introduced animals moved a significant distance away from the introduction site ",
                      "during the establishment phase. The details and destination of post-release dispersal is ",
                      "impossible to predict, but the user can set the probability of otters dispersing."),
               placement = "top", trigger = "hover")
    #
    output$AgeMvAdj_slider <- renderUI({
        sliderInput("AgeMvAdj","Dispersal probability adjustment: subadults relative to adults", 
                    min = 0, max = 1, value = 0.5, step = .01, round = TRUE)
    })
    addPopover(session,"AgeMvAdj_slider",
               "Dispersal probability for subadults relative to adults: more information",
               paste0("In previous sea otter translocations it has been observed that subadult animals ",
                      "may be less likely to disperse (i.e. more likely to remain near the introduction site). ",
                      "This parameter adjusts the likelihood of dispersal for subadults compared to adults: ",
                      "a value of 0.25 would mean that subadults are 1/4 as likely to disperse as adults."),
               placement = "top", trigger = "hover")
    #
    output$EstryMvAdj_slider <- renderUI({
        sliderInput("EstryMvAdj","Dispersal probability adjustment: otters in estuaries", 
                    min = 0, max = 1, value = 0.5, step = .01, round = TRUE)
    })
    addPopover(session,"EstryMvAdj_slider",
               "Dispersal probability for otters in estuaries: more information",
               paste0("Based on several lines of evidence it has been suggested that otters re-introduced to estuaries",
                      "may be less likely to disperse (i.e. more likely to remain near release sites) than otters added ",
                      "to outer coast habitats. This parameter adjusts the likelihood of dispersal for estuaries ",
                      "compared to open coast: a value of 0.25 means otters are 1/4 as likely to disperse post-introduction."),
               placement = "top", trigger = "hover")
    #
    output$Mvmort_slider <- renderUI({
        sliderInput("Mvmort","Mortality rate of otters that disperse during establishment phase ", 
                    min = 0, max = 1, value = 0.75, step = .01, round = TRUE)
    })
    addPopover(session,"Mvmort_slider",
               "Mortality for establishment phase dispersers: more information",
               paste0("The fates of otters that disperse away from a re-introduction sites is hard to determine in most cases: ",
                      "in some reintroductions there has been high levels of mortality for dispersers, in others there is emmigration ",
                      "to a different region altogether. This parameter sets the expected loss-rate for dispersers: the proportion ",
                      "that die or move entirely out of the study area (and are effectively lost to the Oregon meta-population."),
               placement = "top", trigger = "hover")
    #
    output$AddMortE_slider <- renderUI({
        sliderInput("AddMortE","Excess annual mortality during establishment phase", 
                    min = 0, max = .5, value = .14, step = .01, round = TRUE)
    })
    addPopover(session,"AddMortE_slider",
               "Excess mortality during establishment: more information",
               paste0("During the establishment phase of an introduced population, there may be higher than ",
                      "average levels of mortality as the introduced animals become accustomed to their new habitat. ",
                      "In past translocations, excess annual mortality rates of 0.1 - 0.25 have caused translocated ",
                      "populations to decline substantially during the establishment phase"),
               placement = "top", trigger = "hover")
    #
    output$V_sp_slider <- renderUI({
        sliderInput("V_sp","Average expected rate of range exansion (km/yr)", 
                    min = .5, max = 8, value = 2, step = .1, round = FALSE)
    })    
    addPopover(session,"V_sp_slider",
               "Rate of range expansion: more information",
               paste0("Habitat use by the initial sea otter population will likely be limited to a relatively small area of the coast. ",
                      "As this initial population grows, it's distribution (range of occupancy) will spread outwards along the coastline, ",
                      "encompassing more habitat. The rate of range expansion is measured as the speed at which the population front moves ",   
                      "along the coastline. In other populations, this expansion speed has varied from 1-5 km/year."),
               placement = "top", trigger = "hover")
    #
    output$Rmax_slider <- renderUI({
        sliderInput("Rmax","Maximum annual growth rate (at low density)", 
                    min = .1, max = .25, value = .18, step = .01, round = FALSE)
    })      
    addPopover(session,"Rmax_slider",
               "Maximum annual growth rate: more information",
               paste0("Sea otter populations tend to show the highest rate of growth at low densities: as local abundance increases, ", 
                      "the growth rate slows until it eventually reaches 0 when population abundance reaches 'K'. This parameter allows the user ",
                      "to adjust the maximum rate of growth (at low densities): in most sea otter populations this value is close to 0.2"),
               placement = "top", trigger = "hover")    
    #
    output$Estoch_slider <- renderUI({
        sliderInput("Estoch","Environmental stochasticity (SD in annual growth rate)", 
                    min = .01, max = .3, value = .1, step = .01, round = FALSE)
    })  
    addPopover(session,"Estoch_slider",
               "Environmental stochasticity: more information",
               paste0("The average rate of growth for a recovering sea otter population in a given area can be predicted as a function of ",
                      "the local density with respect to carrying capacity, or 'K'. However, year-to-year variation in environmental ",
                      "conditions and prey population dynamics can lead to unpredictable deviations in growth rate, referred to as 'environmental ",
                      "stochasticity'. This parameter controls the degree of variation in growth rates: typical values are 0.05 - 0.15"),
               placement = "top", trigger = "hover")       
    #
    output$Theta_slider <- renderUI({
        sliderInput("Theta","Value of 'theta' (for theta-logistic growth)", 
                    min = .5, max = 2, value = .9, step = .1, round = FALSE)
    })  
    addPopover(session,"Theta_slider",
               "Value of 'theta': more information",
               paste0("The average rate of growth for a recovering sea otter population in a given area can be predicted as a function of ",
                      "the local density with respect to carrying capacity, or 'K'. One of the parameters of this function is 'theta', ",
                      "which determines the onset of reduced growth rates at higher densities: values <1 lead to onset of reduced ",
                      "growth rates at fairly low densities, while values >1 mean that reductions in growth occur only at high densities."),
               placement = "top", trigger = "hover")        
    # 
})
# Create Shiny app ----
shinyApp(ui = ui, server = sv)
