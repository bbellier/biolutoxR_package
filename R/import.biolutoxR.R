# import.biolutoxR() function
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#
# USeful Packages
#
#' @import shiny
#' @import shinyMatrix
#' @import shinythemes
#' @import tidyverse
#' @import drc
#' @import ed50
#' @import openxlsx



# ############################################################################ #
# -------------------------- FONCTION import.biolutoxR() -----------------------
# ############################################################################ #



import.biolutoxR <- function() {
  
  
  
  # ######################################################################## #
  # ----  -1- Chargement des packages  ---------------------------------------
  # ######################################################################## #
  
  
  
  # ---- Liste des packages requis ------------------------------------------- #
  
  
  
  # Liste des packages requis
  required_packages <- c(
    "shiny",
    "shinyMatrix",
    "shinythemes",
    "tidyverse",
    "drc",
    "ed50",
    "openxlsx",
    "remotes"
  )
  # -------------------------------------------------------------------------- #
  
  
  
  # ---- Installation & Chargement des packages ------------------------------ #
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  for (package in required_packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      install.packages(package)
    } else if (package == "ggplot2" && packageVersion("ggplot2") != "3.5.1") {
      remotes::install_version("ggplot2", version = "3.5.1")
    }
    library(package, character.only = TRUE)
  }
  # -------------------------------------------------------------------------- #
  
  
  
  # ---- Suppression des objets inutiles ------------------------------------- #
  rm(package)
  rm(required_packages)
  # -------------------------------------------------------------------------- #
  
  
  
  # ######################################################################## #
  # ----  -2- Application Shiny  ---------------------------------------------
  # ######################################################################## #
  
  
  
  ## A. Style ####
  
  
  
  # ---- Style général de l'application -------------------------------------- #
  style <-
    ".container {
    text-align: center;
  }
  .matrice-container {
    display: inline-block;
  }
  .cellInput {
    width: 205px;
    margin: 5px;
    display: inline-block;
  }
  .shiny-input-container {
    margin: 0px;
  }
  .row-label {
    width: 20px;
    display: inline-block;
    text-align: right;
    margin: 5px;
    font-weight: bold;
  }
  .col-label {
    width: 205px;
    display: inline-block;
    text-align: center;
    margin: 5px;
    font-weight: bold;
    transform: translateX(25px);
  }
  .empty-cell {
    width: 30px;
  }
  .shiny-input-container {
    margin: 0px;
  }"
  
  
  # -------------------------------------------------------------------------- #
  
  
  
  ## B. Fonctions ####
  
  
  
  # ---- Fonction pour calculer une ECx -------------------------------------- #
  toxicity_data_ecX <- function(data, X = 50) {
    drm_model <- drm(mean_perc_inhib_corr ~ as.numeric(dil), data = data,
                     fct = LL.4(fixed=c(NA,0,100,NA),names=c("Slope","Lower Limit","Upper Limit","ED50")))
    ec <- data.frame(ED(drm_model, c(X), interval = "delta"))
    ec_value <- paste0(round(ec$`Estimate`,2), " +/- ", round(ec$Std..Error,2), " (95% confidence interval: ", round(ec[1, "Lower"], 2), "-", round(ec[1, "Upper"], 2), ")")
    assign("ec", ec, envir = .GlobalEnv)
    return(ec_value)
  }
  # -------------------------------------------------------------------------- #
  
  
  
  # ---- Fonction pour crééer le tableau final ------------------------------- #
  toxicity_data_table <- function(data, substance = NA) {
    table_ec <- data %>%
      filter(!is.na(perc_inhib_corr)) %>%
      filter(sol == as.character(substance)) %>%
      group_by(sol, dil) %>%
      summarise(mean_perc_inhib_corr = mean(perc_inhib_corr, na.rm = TRUE)) %>%
      ungroup()
    return(table_ec)
  }
  # -------------------------------------------------------------------------- #
  
  
  
  # ---- Fonction pour représenter graphiquement une ECx --------------------- #
  toxicity_data_table_plot <- function(data) {
    drm_model <- drm(as.numeric(mean_perc_inhib_corr) ~ as.numeric(dil), data = data,
                     fct = LL.4(fixed = c(NA, 0, 100, NA),
                                names = c("Slope", "Lower Limit", "Upper Limit", "EC50")))
    summary(drm_model)
    summary_drm <- summary(drm_model)
    summ_signf <- data.frame(summary_drm$coefficients) %>% rownames_to_column("type") %>%
      mutate(p_value = ifelse(p.value < 0.05, "< 0.05", "> 0.05")) %>%
      mutate(type = ifelse(type == "Slope:(Intercept)", "Slope", "EC50")) %>%
      dplyr::select(c(type, p_value, Estimate)) %>%
      mutate(Estimate = round(Estimate, 2))
    newdata <- expand.grid(conc=exp(seq(log(0.5), log(100), length=100)))
    pm <- predict(drm_model, newdata=newdata, interval="confidence")
    newdata$p <- pm[,1]
    newdata$pmin <- pm[,2]
    newdata$pmax <- pm[,3]
    data$dil0 <- as.numeric(data$dil)
    data$dil0[data$dil0 == 0] <- 0.5
    predicted <- predict(drm_model)
    residuals <- data$mean_perc_inhib_corr - predicted
    SSE <- sum(residuals^2)
    SST <- sum((data$mean_perc_inhib_corr - mean(data$mean_perc_inhib_corr))^2)
    R2 <- 1 - SSE/SST
    coef_model <- coef(drm_model)
    equation <- paste0("y = ", "~0", " + ", "(~100 - ~0)", " / (1 + exp(", round(coef_model[1], 2), " * (log(x) - log(", round(coef_model[2], 2), ")", ")))")
    ggplot(data, aes(x = dil0, y = mean_perc_inhib_corr)) +
      geom_point(shape = 1, size = 2.5) +
      geom_ribbon(data=newdata, aes(x=conc, y=p, ymin=pmin, ymax=pmax), alpha=0.2, fill = "#b1c9ec") +
      geom_line(data=newdata, aes(x=conc, y=p), color = "blue") +
      theme_test() +
      theme(axis.text = element_text(color = "black", size = 12, face = "italic"),
            axis.title = element_text(color = "black", size = 18),
            axis.title.x = element_text(margin = margin(t = 13)),
            axis.title.y = element_text(margin = margin(r = 13)),
            legend.text = element_text(color = "black", size = 12),
            legend.title = element_text(color = "black", size = 12),
            panel.border = element_rect(colour = "black", fill=NA, size=1)) +
      coord_trans(x="log") +
      annotate("text", x = 1, y = 100, label = paste("R² =", round(R2, 2)), size = 5, hjust = 0) +
      annotate("text", x = 1, y = 90, label = "Dose-response curve equation:", size = 4, hjust = 0) +
      annotate("text", x = 1, y = 85, label = equation, size = 4, hjust = 0) +
      annotate("text", x = 1, y = 80, label = paste0("  - ", summ_signf[1,]$type, " (estimation: ", summ_signf[1,]$Estimate, ")", " with a p-value ", summ_signf[1,]$p_value), size = 4, hjust = 0) +
      annotate("text", x = 1, y = 75, label = paste0("  - ", summ_signf[2,]$type, " (estimation: ", summ_signf[2,]$Estimate, ")", " with a p-value ", summ_signf[2,]$p_value), size = 4, hjust = 0) +
      labs(x = "Dilution (in %)", y = "Biolu. inhibition (in %)") +
      scale_y_continuous(limits = c(-10, 110), breaks = c(0, 20, 40, 60, 80, 100))+
      scale_x_continuous(limits = c(1, 100), breaks = c(1.5625, 3.125, 6.25, 12.5, 25, 50, 100))
  }
  # -------------------------------------------------------------------------- #
  
  
  
  # ---- Fonction pour nettoyer l'environnement de travail ------------------- #
  clean_workspace <- function() {
    objects_to_remove <- data.frame(obj = ls())
    for (obj in objects_to_remove) {
      if (exists(obj, envir = .GlobalEnv)) {
        rm(list = obj, envir = .GlobalEnv)
      }
    }
  }
  # -------------------------------------------------------------------------- #
  
  
  
  # ---- Fonction pour filtrer les données "time" ---------------------------- #
  filter_data <- function(data, filter_t1, filter_t2, filter_t3) {
    filtered_data <- data
    if (!filter_t1) {
      filtered_data <- filtered_data %>% filter(time != "Time 1")
    }
    if (!filter_t2) {
      filtered_data <- filtered_data %>% filter(time != "Time 2")
    }
    if (!filter_t3) {
      filtered_data <- filtered_data %>% filter(time != "Time 3")
    }
    return(filtered_data)
  }
  # -------------------------------------------------------------------------- #
  
  
  
  # ---- Fonction pour filtrer les données "dil" ----------------------------- #
  filter_data_dil <- function(data, selected_dilutions) {
    filtered_data <- data %>% filter(dil %in% selected_dilutions)
    return(filtered_data)
  }
  # -------------------------------------------------------------------------- #
  
  
  
  ## C. User Interface (UI) ####
  
  
  
  # ---- UI part ------------------------------------------------------------- #
  
  biolutoxR_ui <-
    
  # --> FluidPage
  fluidPage(
    
    # --> Définition du style (visuel) de l'application
    theme = shinytheme("yeti"),
    tags$head(tags$style(HTML(style))),
    
    # --> Création d'une barre de navigation
    navbarPage("biolutoxR",
               
               # --> Panel "Présentation"
               tabPanel("Presentation",
                        
                        h3(HTML("<b>Welcome to the biolutoxR app!</b>")),
                        h5(HTML("<i>This R-Shiny application facilitates data analysis for toxicity tests based on bacterial bioluminescence inhibition.</i>")),
                        
                        br(),
                        br(),
                        
                        h3(HTML("<b>What is a toxicity test based on bacterial bioluminescence inhibition?</b>")),
                        p("Bacterial bioassays using bioluminescence inhibition are used to assess the bioluminescent response of bacteria to exposure to a solution of interest at different exposure times, typically 5, 15 and 30 minutes. Under optimal conditions, these bacteria produce bioluminescence. This bioluminescence is directly linked to the respiratory metabolic process of the bacteria. However, when bacteria are exposed to a toxic substance, this metabolic process is disrupted, leading to an inhibition of bioluminescence. This reaction is therefore exploited in this type of test to study the direct link between light intensity and the level of toxicity of a sample of interest compared with a control sample."),
                        p("An image for resume the manipulation:"),   
                        
                        tags$a(
                          "Link to the image",
                          target = "_blank",
                          href = "https://github.com/bbellier/biolutoxR_package/blob/main/img/img_white_bckgrnd.png"
                        ),
                        
                        br(),
                        br(),
                        
                        h3(HTML("<b>How use this app?</b>")),
                        h4(HTML("<i>> Settings tab:</i>")),
                        p("- Download your '.xlsx' file (in the correct format)."),
                        p("- Complete the following fields: Manipulation name (name of the experiment), Negative control name (for example NegControl), Time number (for example 5, 15 and 30)."),
                        h4(HTML("<i>> Table tab:</i>")),
                        p("- This tab print the cleaned data."),
                        p("- The cleaned final data table can be exported in '.csv' or '.xlsx' format."),
                        h4(HTML("<i>> Plot tab:</i>")),
                        p("- To visualize test results."),
                        p("- To print dose-response curve."),
                        p("- To calculate ECx values."),
                        h4(HTML("<i>> Exit tab:</i>")),
                        p("- To quit this application."),
                        
                        p("An image for resume the tab functions:"),   
                        
                        tags$a(
                          "Link to the image",
                          target = "_blank",
                          href = "https://github.com/bbellier/biolutoxR_package/blob/main/img/img_pres_tabs_2.png"
                        ),
                        
                        br(),
                        
                        h3(HTML("<b style='color: red;'>Warning</b>")),
                        p(HTML("<b style='color: red;'> - Verify the format of your .xlsx format.</b>")),
                        p(HTML("<b style='color: red;'> - Open the tabs in order to clean the data input, print the plots and obtain the reference ecotoxicity data.</b>")),
                        
                        br(),
                        br(),
                        br(),
                        
                        h5(HTML("<i style='text-align: right; display: block;'>A R-Shiny app package developed by Bellier Benjamin and Le Picard Coralie.</i>"))
                        
               ), # Fermeture du panel "Présentation"
               
               # --> Panel "Settings"
               tabPanel("Settings",
                        
                        h3(HTML("<b>Data settings:</b>")),
                        
                        br(),
                        
                        fileInput("file", "Download your .xlsx file", accept = ".xlsx"),
                        textInput("manip_name", "Manipulation name:"),
                        textInput("neg", "Negative control name", value = "NegControl"),
                        textInput("t_input_t1", "Time 1 (in min):"),
                        textInput("t_input_t2", "Time 2 (in min):"),
                        textInput("t_input_t3", "Time 3 (in min):"),

                        br()
               ), # Fermeture du panel "Settings"
               
               # --> Panel "Table"
               tabPanel("Table",
                        
                        h3(HTML("<b>Final data table:</b>")),
                        downloadButton("download_table_csv", "Download CSV"),
                        downloadButton("download_table_xlsx", "Download XLSX"),
                        
                        br(),
                        br(),
                        
                        dataTableOutput("final_data_table")
                        
               ), # Fermeture du panel "Table"
               
               # --> Panel "Plot"
               tabPanel("Plot",
                        
                        h3(HTML("<b>Final data plot:</b>")),
                        
                        fluidRow(
                          column(
                            width = 4,
                            h4(HTML("Settings:")),
                            h5(HTML("Select time(s) to conserve:")),
                            fluidRow(
                              column(2, align = "left", uiOutput("dynamic_filter_t1")),
                              column(2, align = "left", uiOutput("dynamic_filter_t2")),
                              column(2, align = "left", uiOutput("dynamic_filter_t3"))
                            ),
                            
                            br(),
                            
                            h5(HTML("Select dilution(s) to conserve:")),
                            uiOutput("dynamic_dilution_checkboxes"),
                            
                            br(),
                            
                            h5(HTML("Select variables:")),
                            selectInput("x_variable", "Variable X",
                                        choices = c("sol", "dil", "rep_bio", "sol_sh", "id", "time", "set", "manip_name")),
                            selectInput("y_variable", "Variable Y",
                                        choices = c("biolu", "biolu_mean_neg_control", "biolu_corr", "perc_inhib", "perc_inhib_corr")),
                            selectInput("color_variable", "Color",
                                        choices = c("sol", "dil", "rep_bio", "sol_sh", "id", "time", "set", "manip_name"))
                          ),
                          column(
                            width = 8,
                            plotOutput("plot_biolutoxR"),
                            
                            br(),
                            
                            fluidRow(
                              column(3),
                              column(6, div(style = "text-align: center;", uiOutput("dynamic_substance_checkboxes"))),
                              column(3)
                            ),
                            
                            br()
                            
                          )
                        ),
                        
                        br(),
                        
                        h3(HTML("<b>Dose-Response curve:</b>")),
                        fluidRow(
                          column(
                            width = 4,
                            h4(HTML("Settings:")),
                            h5(HTML("Select time(s) to conserve:")),
                            fluidRow(
                              column(2, align = "left", uiOutput("dynamic_filter_t1_2")),
                              column(2, align = "left", uiOutput("dynamic_filter_t2_2")),
                              column(2, align = "left", uiOutput("dynamic_filter_t3_2"))
                            ),
                            
                            br(),
                            
                            h5(HTML("Select the variable of interest:")),
                            selectInput("sol", "Variable:", choices = NULL, selectize = FALSE),
                            
                            br(),
                            
                            h5(HTML("Select the % of efectiveness")),
                            numericInput("ecX_input_value", "Value (in %)", 50, min = 1, max = 100),
                            
                            br(),
                            
                            h5(HTML("Result:")),
                            textOutput("ecX_text"),
                            
                            br(),
                            
                            h5(HTML("Explanation:")),
                            textOutput("ecX_text2"),
                            
                            br()
                            
                          ),
                          column(
                            width = 8,
                            plotOutput("ec50_plot"),
                            h6(HTML("The package automatically calculates EC50 values with at least one value per dilution, but please ensure that your data are relevant (e.g. number of dilutions, biological replicates and set)"), style = "text-align: center; color: red; font-weight: bold;")
                          )
                        )
               ), # Fermeture du panel "Plot"
               
               # --> Panel "Exit"
               tabPanel("Exit",
                        actionButton("exit_btn", "Exit")
               ) # Fermeture du panel "Exit"
               
    ) # Fermeture du style
    
  ) # Fermeture de la FluidPage
  # -------------------------------------------------------------------------- #
  
  
  
  ## D. Server ####
  
  
  
  # ---- Server part --------------------------------------------------------- #
  biolutoxR_server <- function(input, output, session) {
    
    
    
  ### a. Définition des valeurs par défaut ####
    
    
    
  # Données traitées  
  data <- reactive({
    
    req(input$file)
    
    scheme <- read.xlsx(input$file$datapath, sheet = "scheme")
    time1 <- read.xlsx(input$file$datapath, sheet = "time 1")
    time2 <- read.xlsx(input$file$datapath, sheet = "time 2")
    time3 <- read.xlsx(input$file$datapath, sheet = "time 3")
    
    neg <- input$neg
    
    combined_data <- bind_rows(
      scheme %>% left_join(time1, by = c("set", "row", "column")) %>% mutate(time = "Time 1") %>% mutate(biolu = as.numeric(biolu)),
      scheme %>% left_join(time2, by = c("set", "row", "column")) %>% mutate(time = "Time 2") %>% mutate(biolu = as.numeric(biolu)),
      scheme %>% left_join(time3, by = c("set", "row", "column")) %>% mutate(time = "Time 3") %>% mutate(biolu = as.numeric(biolu))
    ) %>%
      mutate(biolu = as.numeric(biolu)) %>%
      group_by(time, set) %>%
      mutate(
        biolu_mean_neg_control = round(mean(biolu[sol == neg], na.rm = TRUE), 2),
        biolu_corr = round(biolu * 100 / biolu_mean_neg_control, 2),
        perc_inhib = round(100 - biolu_corr, 2),
        perc_inhib_corr = ifelse(perc_inhib > 100, 100, ifelse(perc_inhib < 0, 0, round(perc_inhib, 2)))
      ) %>% 
      mutate(perc_inhib_corr = as.numeric(perc_inhib_corr)) %>% 
      mutate(id = paste0(row, column)) %>% 
      mutate(set = round(set, 0)) %>% 
      mutate(manip_name = input$manip_name) %>% 
      dplyr::select(manip_name, set, id, sol, dil,	rep_bio,	sol_sh, time,
                    biolu, biolu_mean_neg_control, biolu_corr, perc_inhib, perc_inhib_corr)
    combined_data
    
  })
  
  
  
  # Valeur de "manip_name" pour l'exemple
  observe({
    if (is.null(input$manip_name) || input$manip_name == "") {
      updateTextInput(session, "manip_name", value = "MyManip")
    }
  })
  
  
  
  # Valeur de T1 pour l'exemple
  observe({
    if (is.null(input$t_input_t1) || input$t_input_t1 == "") {
      updateTextInput(session, "t_input_t1", value = "5")
    }
  })
  output$dynamic_filter_t1 <- renderUI({
    checkboxInput("filter_t1", input$t_input_t1, value = FALSE)
  })
  output$dynamic_filter_t1_2 <- renderUI({
    checkboxInput("filter_t1_2", input$t_input_t1, value = FALSE)
  })
  
  
  
  # Valeur de T2 pour l'exemple
  observe({
    if (is.null(input$t_input_t2) || input$t_input_t2 == "") {
      updateTextInput(session, "t_input_t2", value = "15")
    }
  })
  output$dynamic_filter_t2 <- renderUI({
    checkboxInput("filter_t2", input$t_input_t2, value = FALSE)
  })
  output$dynamic_filter_t2_2 <- renderUI({
    checkboxInput("filter_t2_2", input$t_input_t2, value = FALSE)
  })
  
  
  
  # Valeur de T3 pour l'exemple
  observe({
    if (is.null(input$t_input_t3) || input$t_input_t3 == "") {
      updateTextInput(session, "t_input_t3", value = "30")
    }
  })
  output$dynamic_filter_t3 <- renderUI({
    checkboxInput("filter_t3", input$t_input_t3, value = TRUE)
  })
  output$dynamic_filter_t3_2 <- renderUI({
    checkboxInput("filter_t3_2", input$t_input_t3, value = TRUE)
  })
  
  
  
  # Valeur de "neg_control" pour l'exemple
  observe({
    if (is.null(input$neg_control) || input$neg_control == "") {
      updateTextInput(session, "neg_control", value = "NegControl")
    }
  })
  
  
  
  # Valeur de "sol" pour l'exemple
  observe({
    if (is.null(input$sol) || input$sol == "") {
      updateTextInput(session, "sol", value = "Substance A")
    }
  })
  
  
  
  ### b. Création d'un tableau final ####
  
  
  
  # Réalisation d'un tableau final
  output$final_data_table <-
    
    renderDataTable(
      
      {
        
        final_data_corrected <- data()
        final_data_corrected <- final_data_corrected %>% filter(!is.na(biolu))
        
        assign("final_data_corrected", final_data_corrected, envir = .GlobalEnv)
        
        return(final_data_corrected)
        
      },
      
      options = list(
        lengthMenu = c(30, 60, 90),
        pageLength = 30,
        scrollY = "300px",
        scrollX = TRUE
        
      )
      
    )
  
  
  
  # Exportation de la table en formlat ".csv"
  output$download_table_csv <- downloadHandler(
    filename = function() {"final_data_table.csv"},
    content = function(file) {write.csv(data(), file, row.names = FALSE)}
  )
  
  
  
  # Exportation de la table en format ".xlsx"
  output$download_table_xlsx <- downloadHandler(
    filename = function() {"final_data_table.xlsx"},
    content = function(file) {
      openxlsx::write.xlsx(data(), file, rowNames = FALSE)
    }
  )
  
  
  
  ### c. Création d'un graphique final "général" ####
  
  
  
  # Adapter les données pour le graphique final
  final_data_corrected <-
    
    reactive(
      
      {
        
        observe({
          data <- final_data_corrected() %>% 
            filter(sol != input$neg)
          sol_names <- unique(data$sol)
          updateSelectInput(session, "sol", choices = sol_names)
        })
        
        
        final_data_corrected <- data()
        final_data_corrected <- final_data_corrected %>% filter(!is.na(biolu))
        final_data_corrected <- final_data_corrected %>% arrange(dil)
        
        niveaux <- sort(as.numeric(unique(final_data_corrected$dil)))
        
        return(final_data_corrected)
        
      })
  
  
  
  # Sélectionner certaines variables par défaut pour le graphique
  observe({
    updateSelectInput(session, "x_variable", choices = c("sol", "dil", "rep_bio", "sol_sh", "id", "time", "set", "manip_name"), selected = "sol")
    updateSelectInput(session, "y_variable", choices = c("biolu", "biolu_mean_neg_control", "biolu_corr", "perc_inhib", "perc_inhib_corr"), selected = "perc_inhib_corr")
    updateSelectInput(session, "color_variable", choices = c("sol", "dil", "rep_bio", "sol_sh", "id", "time", "set", "manip_name"), selected = "dil")
  })
  
  
  
  # Choix des substances dynamique en fonction des valeurs entrées par l'utilisateur
  output$dynamic_substance_checkboxes <-
    
    renderUI({
      
      unique_substances <- unique(final_data_corrected()$sol)
      
      checkbox_group <- checkboxGroupInput(
        inputId = "selected_substances",
        label = "Select substances to focus on:",
        choices = unique_substances,
        selected = unique_substances,
        inline = TRUE
      )
      
      return(checkbox_group)
      
    })
  
  
  
  # Choix des dilutions dynamique en fonction des valeurs entrées par l'utilisateur
  output$dynamic_dilution_checkboxes <-
    
    renderUI({
      
      unique_dilutions <- sort(as.numeric(unique(final_data_corrected()$dil)))
      
      checkbox_group <- checkboxGroupInput(
        inputId = "selected_dilutions",
        label = "Dilution(s) to include:",
        choices = unique_dilutions,
        selected = unique_dilutions
      )
      
      return(checkbox_group)
      
    })
  
  
  
  # Adapter les données en fonction des choix pour le graphique
  correct_data_biolutox_filtered <-
    
    reactive({
      
      final_data_corrected <- final_data_corrected()
      final_data_corrected <- filter_data(final_data_corrected, input$filter_t1, input$filter_t2, input$filter_t3)
      selected_dilutions <- as.numeric(input$selected_dilutions)
      final_data_corrected <- filter_data_dil(final_data_corrected, selected_dilutions)
      
      return(final_data_corrected)
      
    })
  
  
  
  # Création du graphique final
  output$plot_biolutoxR <-
    
    renderPlot({
      
      req(input$x_variable, input$y_variable, input$color_variable)
      
      data_plot <- correct_data_biolutox_filtered()
      data_plot$dil = factor(data_plot$dil, levels = sort(as.numeric(unique(data_plot$dil))))
      data_plot$time = factor(data_plot$time, levels = c("Time 1", "Time 2", "Time 3"))
      
      selected_substances <- input$selected_substances
      
      if (!is.null(selected_substances) && length(selected_substances) > 0) {
        data_plot <- data_plot %>% filter(sol %in% selected_substances)
      }
      
      p_data <- ggplot(data_plot, aes_string(x = input$x_variable,
                                             y = input$y_variable,
                                             color = input$color_variable)) +
        geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.85)) +
        geom_jitter(alpha = 0.5, position = position_dodge(width = 0.85)) +
        stat_summary(fun = "mean", geom = "point", shape = 4, size = 3, position = position_dodge(width = 0.85)) +
        theme_test() +
        theme(axis.text = element_text(color = "black", size = 12, face = "italic"),
              axis.title = element_text(color = "black", size = 18),
              axis.title.x = element_text(margin = margin(t = 13)),
              axis.title.y = element_text(margin = margin(r = 13)),
              legend.text = element_text(color = "black", size = 12),
              legend.title = element_text(color = "black", size = 12),
              legend.position = "top",
              legend.justification.top = "right",
              panel.border = element_rect(colour = "black", fill=NA, size=1))
      print(p_data)
      
    })
  
  
  
  ### d. Création d'un graphique final "toxicité" ####
  
  
  
  # Création du graphique final
  output$ec50_plot <-
    
    renderPlot({
      
      tryCatch({
        ec_data <- final_data_corrected()
        ec_data <- filter_data(ec_data, input$filter_t1_2, input$filter_t2_2, input$filter_t3_2)
        table_ec <- toxicity_data_table(ec_data, substance = input$sol)
        toxicity_data_table_plot(table_ec)
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, "NA", cex = 2, col = "red")
      })
      
    })
  
  
  
  # Création du résultat
  output$ecX_text <- renderText({
    
    tryCatch({
      ec_data <- final_data_corrected()
      ec_data <- filter_data(ec_data, input$filter_t1_2, input$filter_t2_2, input$filter_t3_2)
      table_ec <- toxicity_data_table(ec_data, substance = input$sol)
      ec_X_value <- toxicity_data_ecX(table_ec, X = input$ecX_input_value)
      paste0("EC", input$ecX_input_value, " value: ", ec_X_value)
    }, error = function(e) {
      "Changer de substance"
    })
    
  })
  
  
  
  # Création d'une phrase de conclusion
  output$ecX_text2 <- renderText({
    
    tryCatch({
      ec_data <- final_data_corrected()
      ec_data <- filter_data(ec_data, input$filter_t1_2, input$filter_t2_2, input$filter_t3_2)
      table_ec <- toxicity_data_table(ec_data, substance = input$sol)
      ec_X_value <- toxicity_data_ecX(table_ec, X = input$ecX_input_value)
      paste0("Concretely, bacterial bioluminescence production was inhibited by ",  input$ecX_input_value, "% at a dilution percentage of ", ec_X_value, " of the solution of interest (here, ", input$sol, ") after a considered exposure time (depending on the time(s) selected).")
    }, error = function(e) {
      "Changer de substance"
    })
    
  })
  

  
  ## e. Création d'un graphique final "toxicité" ####
  observeEvent(input$exit_btn,{
    stopApp()
    clean_workspace()
  })
  
  
  
  } # Fermeture du server
  # -------------------------------------------------------------------------- #
  
  
  
  ## E. Terminaison de l'application ####
  shinyApp(ui = biolutoxR_ui, server = biolutoxR_server)
  
  
  
} # Fermeture de l'application



# ############################################################################ #
# ------------------------ LANCEMENT import.biolutoxR() -----------------------
# ############################################################################ #



# Run la function import.biolutoxR()
import.biolutoxR()


