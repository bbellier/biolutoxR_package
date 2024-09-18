# run.biolutoxR() function
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
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


# ############################################################################ #
# ------------------------ FONCTION example.biolutoxR() ------------------------
# ############################################################################ #


example.biolutoxR <- function() {



  # ######################################################################## #
  # ----  -1- Chargement des packages  ---------------------------------------
  # ######################################################################## #



  # ---- Liste des packages requis ------------------------------------------- #
  required_packages <- c(
    "shiny",
    "shinyMatrix",
    "shinythemes",
    "tidyverse",
    "drc",
    "ed50"
  )
  # -------------------------------------------------------------------------- #




  # ---- Installation & Chargement des packages ------------------------------ #
  for (package in required_packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      install.packages(package)
    } else {
      library(package, character.only = TRUE)
    }
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



  # ---- Fonction pour générer dynamiquement les entrées de la matrice ------- #
  create_matrix <- function(matrix_prefix, ncol, nrow, input, label_box, default_values_list) {
    matrices_ui <- lapply(1:input$n_set, function(i) {
      matrix_name <- paste0(matrix_prefix, "_", i)
      default_values <- default_values_list[[i]]
      matrix_ui <- tagList(
        h4(matrix_name),
        fluidRow(
          div(class = "empty-cell"),
          lapply(1:ncol, function(col) {
            div(class = "col-label", col)
          })
        ),
        lapply(LETTERS[1:nrow], function(row) {
          fluidRow(
            div(class = "row-label", row),
            lapply(1:ncol, function(col) {
              cell_id <- paste(matrix_name, row, col, sep = "_")
              default_val <- default_values[[paste0(row, col)]]
              div(class = "cellInput", textInput(inputId = cell_id, label = NULL, placeholder = label_box, value = default_val))
            })
          )
        })
      )
      div(class = "matrice-container", matrix_ui)
    })
    tagList(matrices_ui)
  }
  # -------------------------------------------------------------------------- #



  # ---- Fonction pour corriger les données finales -------------------------- #
  correct_final_data <- function(data, neg) {
    correct_data_microtox <- data %>%
      mutate(biolu = as.numeric(biolu)) %>%
      left_join(
        data %>%
          mutate(biolu = as.numeric(biolu)) %>%
          filter(sol_type == as.character(neg)) %>%
          group_by(time, set) %>%
          summarise(mean_biolu_nacl = round(mean(biolu, na.rm = TRUE), 2), .groups = "drop"),
        by = c("time", "set")) %>%
      mutate(biolu_corr = biolu*100/mean_biolu_nacl,
             biolu_corr = round(biolu_corr, 2),
             perc_inhib = 100-biolu_corr,
             perc_inhib = round(perc_inhib, 2)) %>%
      mutate(perc_inhib_corr = ifelse(perc_inhib > 100, 100, perc_inhib),
             perc_inhib_corr = ifelse(perc_inhib < 0, 0, perc_inhib),
             perc_inhib_corr = round(perc_inhib_corr, 2))
    assign("correct_data_microtox", correct_data_microtox, envir = .GlobalEnv)
  }
  # -------------------------------------------------------------------------- #



  # ---- Fonction pour crééer le tableau final ------------------------------- #
  toxicity_data_table <- function(data, substance = NA) {
    table_ec <- data %>%
      filter(!is.na(perc_inhib_corr)) %>%
      filter(sol_type == as.character(substance)) %>%
      group_by(sol_type, dil) %>%
      summarise(mean_perc_inhib_corr = mean(perc_inhib_corr, na.rm = TRUE)) %>%
      ungroup()
    return(table_ec)
  }
  # -------------------------------------------------------------------------- #



  # ---- Fonction pour calculer une ECx -------------------------------------- #
  toxicity_data_ecX <- function(data, X = 50) {
    drm_model <- drm(mean_perc_inhib_corr ~ as.numeric(dil), data = data,
                     fct = LL.4(fixed=c(NA,0,100,NA),
                                names=c("Slope","Lower Limit","Upper Limit","ED50")))
    ec <- data.frame(ED(drm_model, c(X)))
    ec_value <- paste0(round(ec$`Estimate`,2), " +/- ", round(ec$Std..Error,2))
    assign("ec", ec, envir = .GlobalEnv)
    return(ec_value)
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
      coord_trans(x="log") +
      annotate("text", x = 1, y = 105, label = paste("R² =", round(R2, 2)), size = 5, hjust = 0) +
      annotate("text", x = 1, y = 95, label = "Equation:", size = 4, hjust = 0) +
      annotate("text", x = 1, y = 90, label = equation, size = 4, hjust = 0) +
      annotate("text", x = 1, y = 85, label = paste0("  - ", summ_signf[1,]$type, " (Estimate:", summ_signf[1,]$Estimate, ")", " with a p-value ", summ_signf[1,]$p_value), size = 4, hjust = 0) +
      annotate("text", x = 1, y = 80, label = paste0("  - ", summ_signf[2,]$type, " (Estimate:", summ_signf[2,]$Estimate, ")", " with a p-value ", summ_signf[2,]$p_value), size = 4, hjust = 0) +
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
  filter_data <- function(data, filter_t5, filter_t15, filter_t30) {
    filtered_data <- data
    if (!filter_t5) {
      filtered_data <- filtered_data %>% filter(time != "Time 5")
    }
    if (!filter_t15) {
      filtered_data <- filtered_data %>% filter(time != "Time 15")
    }
    if (!filter_t30) {
      filtered_data <- filtered_data %>% filter(time != "Time 30")
    }
    return(filtered_data)
  }
  # -------------------------------------------------------------------------- #



  # ---- Fonction pour filtrer les données "dil" ---------------------------- #
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

                          fluidRow(
                            column(
                              width = 1,
                              img(src = "https://cdn.icon-icons.com/icons2/881/PNG/512/Fish_Food_icon-icons.com_68747.png",
                                  width = "80%",
                                  height = "80%")
                            ),
                            column(
                              width = 11,
                              h3(HTML("<b>Welcome to the biolutoxR app!</b>")),
                              p("This app is very useful to facilitate data analysis of toxicity test based on bacterial bioluminescence inhibition data.")
                            )
                          ),

                          br(),

                          h3(HTML("<b>What is a toxicity test based on bacterial bioluminescence inhibition?</b>")),
                          p("Among standardized toxicity tests, toxicity test based on bacterial bioluminescence inhibitions offer the possibility of studying chemical toxicity with 'speed, simplicity, reproducibility, precision, sensitivity, standardization, cost-effectiveness and convenience' (Isenberg, 1993). More concretely, the toxicity test based on bacterial bioluminescence inhibition is a standardized acute toxicity test (bacterial luminescence test ISO 11348-3 and ISO 21338) which evaluates the bioluminescent response of a bacterium to a potentially toxic stress (Johnson, 2005). This bioassay uses Aliivibrio fisheri (formerly Vibrio fisheri), a non-pathogenic gram-negative marine bacterium which, under optimal conditions, produces bioluminescence, a phenomenon directly linked to its respiratory metabolic process. When exposed to a toxic substance, this metabolic process is disrupted, leading to a reduction in bioluminescence. This reaction is therefore exploited in the toxicity test based on bacterial bioluminescence inhibition to study the direct link between light intensity and the toxicity level of a sample compared to a control sample (with no toxic substance)."),

                          br(),

                          h3(HTML("<b>How use this app?</b>")),
                          p("1) Scheme panel: "),
                          p("- Complete the following fields: Manipulation name (name of the experiment), Number of matrices (depending of the set number), Number of lines (often 6), Number of columns (often 5) and Negative control name (for example NaCl)"),
                          p("- Complete matrices, i.e. for each cell, the name of the solution, the dilution, the biological replicate and finally the short name of the solution, all separated by a comma."),
                          p("2) Time X panel: for each time panel, complete matrix based on the results obtained during the experiment."),
                          p("3) Table panel: raw resutls."),
                          p("4) Plot panel: comparison plot and dose-response curve and EC50 and EC20 value."),
                          p("5) Exit panel: to quit this application."),

                          br(),

                          h3(HTML("<b style='color: red;'>Warning</b>")),
                          p(HTML("<b style='color: red;'> - Open the tabs in sequence to load the table and graphics.</b>")),
                          p(HTML("<b style='color: red;'> - Don't add new matrices after you have filled your matrices.</b>"))

                 ), # Fermeture du panel "Présentation"

                 # --> Panel "Scheme"
                 tabPanel("Scheme",

                          h3(HTML("<b>Data settings:</b>")),

                          br(),

                          fluidRow(
                            column(2, align = "center", textInput("manip_name", "Manipulation name:")),
                            column(2, align = "center", numericInput("n_set", "Number of matrices:", min = 1, max = 20, value = 2)),
                            column(2, align = "center", numericInput("n_row", "Number of lines:", min = 1, max = 20, value = 6)),
                            column(2, align = "center", numericInput("n_col", "Number of columns:", min = 1, max = 20, value = 5)),
                            column(2, align = "center", textInput("neg_control", "Negative control name:"))
                          ),
                          h3(HTML("<b>Scheme matrix/matrices:</b>")),
                          uiOutput("matrix_scheme")

                 ), # Fermeture du panel "Scheme"

                 # --> Panel "Time 5"
                 tabPanel("Time 5",

                          h3(HTML("<b>Time 5 matrix/matrices:</b>")),
                          uiOutput("matrix_t5")

                 ), # Fermeture du panel "Time 5"

                 # --> Panel "Time 15"
                 tabPanel("Time 15",

                          h3(HTML("<b>Time 15 matrix/matrices:</b>")),
                          uiOutput("matrix_t15")

                 ), # Fermeture du panel "Time 15"

                 # --> Panel "Time 30"
                 tabPanel("Time 30",

                          h3(HTML("<b>Time 30 matrix/matrices:</b>")),
                          uiOutput("matrix_t30")

                 ), # Fermeture du panel "Time 30"

                 # --> Panel "Table"
                 tabPanel("Table",

                          h3(HTML("<b>Final data table:</b>")),
                          downloadButton("download_table_csv", "Download CSV"),

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
                              h5(HTML("Select time to conserve:")),
                              fluidRow(
                                column(2, align = "left", checkboxInput("filter_t5", "T5", value = FALSE)),
                                column(2, align = "left", checkboxInput("filter_t15", "T15", value = FALSE)),
                                column(2, align = "left", checkboxInput("filter_t30", "T30", value = TRUE))
                              ),

                              br(),

                              h5(HTML("Select dilution to conserve:")),
                              uiOutput("dynamic_dilution_checkboxes"),

                              br(),

                              h5(HTML("Select variables:")),
                              selectInput("x_variable", "Variable X",
                                          choices = c("sol_type", "dil", "rep_bio", "name", "id", "time", "set", "manip_name")),
                              selectInput("y_variable", "Variable Y",
                                          choices = c("biolu", "mean_biolu_nacl", "biolu_corr", "perc_inhib", "perc_inhib_corr")),
                              selectInput("color_variable", "Color",
                                          choices = c("sol_type", "dil", "rep_bio", "name", "id", "time", "set", "manip_name"))
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
                              h5(HTML("Select time variables to conserve:")),
                              fluidRow(
                                column(2, align = "left", checkboxInput("filter_t5_2", "T5", value = TRUE)),
                                column(2, align = "left", checkboxInput("filter_t15_2", "T15", value = TRUE)),
                                column(2, align = "left", checkboxInput("filter_t30_2", "T30", value = TRUE))
                              ),

                              br(),

                              h5(HTML("Select the variable of interest:")),
                              selectInput("sol_type", "Variable:", choices = NULL, selectize = FALSE),

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
                              h5(HTML("Are represented in red the EC50 and in yellow the EC20."), style = "text-align: center;"),
                              h6(HTML("The package automatically calculates EC50 values with at least one value per dilution, but please ensure that your data are relevant (e.g. number of replicates and dilutions)"), style = "text-align: center; color: red; font-weight: bold;")
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



    # Valeurs par défaut pour l'onglet "Scheme" pour l'exemple
    default_values_list_scheme <-
      list(
        list(
          "A1" = "NegControl,100,1,Neg",
          "A2" = "NegControl,100,1,Neg",
          "A3" = "NegControl,100,1,Neg",
          "A4" = "NegControl,100,1,Neg",
          "A5" = "NegControl,100,1,Neg",
          "B1" = "Substance A,100,1,A",
          "B2" = "Substance A,50,1,A",
          "B3" = "Substance A,25,1,A",
          "B4" = "Substance A,12.5,1,A",
          "B5" = "Substance A,6.225,1,A",
          "C1" = "Substance B,100,1,B",
          "C2" = "Substance B,50,1,B",
          "C3" = "Substance B,25,1,B",
          "C4" = "Substance B,12.5,1,B",
          "C5" = "Substance B,6.225,1,B",
          "D1" = "NA",
          "D2" = "NA",
          "D3" = "NA",
          "D4" = "NA",
          "D5" = "NA",
          "E1" = "NA",
          "E2" = "NA",
          "E3" = "NA",
          "E4" = "NA",
          "E5" = "NA",
          "F1" = "NA",
          "F2" = "NA",
          "F3" = "NA",
          "F4" = "NA",
          "F5" = "NA"
        ),
        list(
          "A1" = "NegControl,100,2,Neg",
          "A2" = "NegControl,100,2,Neg",
          "A3" = "NegControl,100,2,Neg",
          "A4" = "NegControl,100,2,Neg",
          "A5" = "NegControl,100,2,Neg",
          "B1" = "Substance A,100,2,A",
          "B2" = "Substance A,50,2,A",
          "B3" = "Substance A,25,2,A",
          "B4" = "Substance A,12.5,2,A",
          "B5" = "Substance A,6.225,2,A",
          "C1" = "Substance B,100,2,B",
          "C2" = "Substance B,50,2,B",
          "C3" = "Substance B,25,2,B",
          "C4" = "Substance B,12.5,2,B",
          "C5" = "Substance B,6.225,2,B",
          "D1" = "NA",
          "D2" = "NA",
          "D3" = "NA",
          "D4" = "NA",
          "D5" = "NA",
          "E1" = "NA",
          "E2" = "NA",
          "E3" = "NA",
          "E4" = "NA",
          "E5" = "NA",
          "F1" = "NA",
          "F2" = "NA",
          "F3" = "NA",
          "F4" = "NA",
          "F5" = "NA"
        )
      )



    # Valeurs par défaut pour l'onglet "Time 5" pour l'exemple
    default_values_list_t5 <-
      list(
        list(
          "A1" = "100",
          "A2" = "96",
          "A3" = "100",
          "A4" = "95",
          "A5" = "90",
          "B1" = "5",
          "B2" = "26",
          "B3" = "32",
          "B4" = "60",
          "B5" = "82",
          "C1" = "2",
          "C2" = "3",
          "C3" = "13",
          "C4" = "26",
          "C5" = "45",
          "D1" = "NA",
          "D2" = "NA",
          "D3" = "NA",
          "D4" = "NA",
          "D5" = "NA",
          "E1" = "NA",
          "E2" = "NA",
          "E3" = "NA",
          "E4" = "NA",
          "E5" = "NA",
          "F1" = "NA",
          "F2" = "NA",
          "F3" = "NA",
          "F4" = "NA",
          "F5" = "NA"
        ),
        list(
          "A1" = "100",
          "A2" = "96",
          "A3" = "100",
          "A4" = "95",
          "A5" = "90",
          "B1" = "5",
          "B2" = "28",
          "B3" = "33",
          "B4" = "64",
          "B5" = "84",
          "C1" = "5",
          "C2" = "6",
          "C3" = "14",
          "C4" = "28",
          "C5" = "50",
          "D1" = "NA",
          "D2" = "NA",
          "D3" = "NA",
          "D4" = "NA",
          "D5" = "NA",
          "E1" = "NA",
          "E2" = "NA",
          "E3" = "NA",
          "E4" = "NA",
          "E5" = "NA",
          "F1" = "NA",
          "F2" = "NA",
          "F3" = "NA",
          "F4" = "NA",
          "F5" = "NA"
        )
      )



    # Valeurs par défaut pour l'onglet "Time 15" pour l'exemple
    default_values_list_t15 <-
      list(
        list(
          "A1" = "100",
          "A2" = "96",
          "A3" = "100",
          "A4" = "95",
          "A5" = "90",
          "B1" = "9",
          "B2" = "31",
          "B3" = "37",
          "B4" = "66",
          "B5" = "87",
          "C1" = "3",
          "C2" = "5",
          "C3" = "23",
          "C4" = "40",
          "C5" = "58",
          "D1" = "NA",
          "D2" = "NA",
          "D3" = "NA",
          "D4" = "NA",
          "D5" = "NA",
          "E1" = "NA",
          "E2" = "NA",
          "E3" = "NA",
          "E4" = "NA",
          "E5" = "NA",
          "F1" = "NA",
          "F2" = "NA",
          "F3" = "NA",
          "F4" = "NA",
          "F5" = "NA"
        ),
        list(
          "A1" = "100",
          "A2" = "96",
          "A3" = "100",
          "A4" = "95",
          "A5" = "90",
          "B1" = "10",
          "B2" = "33",
          "B3" = "39",
          "B4" = "68",
          "B5" = "86",
          "C1" = "4",
          "C2" = "6",
          "C3" = "26",
          "C4" = "43",
          "C5" = "55",
          "D1" = "NA",
          "D2" = "NA",
          "D3" = "NA",
          "D4" = "NA",
          "D5" = "NA",
          "E1" = "NA",
          "E2" = "NA",
          "E3" = "NA",
          "E4" = "NA",
          "E5" = "NA",
          "F1" = "NA",
          "F2" = "NA",
          "F3" = "NA",
          "F4" = "NA",
          "F5" = "NA"
        )
      )



    # Valeurs par défaut pour l'onglet "Time 30" pour l'exemple
    default_values_list_t30 <-
      list(
        list(
          "A1" = "100",
          "A2" = "96",
          "A3" = "100",
          "A4" = "95",
          "A5" = "90",
          "B1" = "15",
          "B2" = "36",
          "B3" = "42",
          "B4" = "69",
          "B5" = "92",
          "C1" = "7",
          "C2" = "9",
          "C3" = "28",
          "C4" = "48",
          "C5" = "64",
          "D1" = "NA",
          "D2" = "NA",
          "D3" = "NA",
          "D4" = "NA",
          "D5" = "NA",
          "E1" = "NA",
          "E2" = "NA",
          "E3" = "NA",
          "E4" = "NA",
          "E5" = "NA",
          "F1" = "NA",
          "F2" = "NA",
          "F3" = "NA",
          "F4" = "NA",
          "F5" = "NA"
        ),
        list(
          "A1" = "100",
          "A2" = "96",
          "A3" = "100",
          "A4" = "95",
          "A5" = "90",
          "B1" = "18",
          "B2" = "35",
          "B3" = "44",
          "B4" = "68",
          "B5" = "94",
          "C1" = "8",
          "C2" = "10",
          "C3" = "31",
          "C4" = "52",
          "C5" = "67",
          "D1" = "NA",
          "D2" = "NA",
          "D3" = "NA",
          "D4" = "NA",
          "D5" = "NA",
          "E1" = "NA",
          "E2" = "NA",
          "E3" = "NA",
          "E4" = "NA",
          "E5" = "NA",
          "F1" = "NA",
          "F2" = "NA",
          "F3" = "NA",
          "F4" = "NA",
          "F5" = "NA"
        )
      )



    # Valeur de "manip_name" pour l'exemple
    observe({
      if (is.null(input$manip_name) || input$manip_name == "") {
        updateTextInput(session, "manip_name", value = "Example")
      }
    })



    # Valeur de "neg_control" pour l'exemple
    observe({
      if (is.null(input$neg_control) || input$neg_control == "") {
        updateTextInput(session, "neg_control", value = "NegControl")
      }
    })



    # Valeur de "sol_type" pour l'exemple
    observe({
      if (is.null(input$sol_type) || input$sol_type == "") {
        updateTextInput(session, "sol_type", value = "Substance A")
      }
    })



    ### b. Modification des matrices ####



    # Adaptation des matrices "Scheme"
    output$matrix_scheme <- renderUI({
      create_matrix("Scheme_Set", input$n_col, input$n_row, input, "sol_type,dil,rep_bio,name", default_values_list_scheme)
    })



    # Adaptation des matrices "Time 5"
    output$matrix_t5 <- renderUI({
      create_matrix("Time 5_Set", input$n_col, input$n_row, input, "biolu", default_values_list_t5)
    })



    # Adaptation des matrices "Time 15"
    output$matrix_t15 <- renderUI({
      create_matrix("Time 15_Set", input$n_col, input$n_row, input, "biolu", default_values_list_t15)
    })



    # Adaptation des matrices "Time 30"
    output$matrix_t30 <- renderUI({
      create_matrix("Time 30_Set", input$n_col, input$n_row, input, "biolu", default_values_list_t30)
    })



    ### c. Extraction des valeurs depuis les matrices ####



    # Extraire les valeurs des matrices
    extract_matrix_values <-
      function(matrix_prefix, ncol, nrow, input) {
        matrix_data <- data.frame(
          matrix = character(),
          row = character(),
          column = integer(),
          value = character(),
          stringsAsFactors = FALSE
          )

        for (i in 1:input$n_set) {
          matrix_name <- paste0(matrix_prefix, "_", i)
          for (row in LETTERS[1:nrow]) {
            for (col in 1:ncol) {
              cell_id <- paste(matrix_name, row, col, sep = "_")
              value <- input[[cell_id]]
              matrix_data <- rbind(matrix_data, data.frame(matrix = matrix_name,
                                                           row = row,
                                                           column = col,
                                                           value = value,
                                                           stringsAsFactors = FALSE))

            }
          }
        }

        return(matrix_data)

      }



    ### d. Création d'un tableau final ####



    # Réalisation d'un tableau final
    output$final_data_table <-

      renderDataTable(

        {

        scheme <- extract_matrix_values("Scheme_Set", input$n_col, input$n_row, input) %>%
          mutate(id = paste0(row, column)) %>%
          rename(scheme_value = value) %>%
          separate(col = matrix, into = c("matrix", "set", "set_value"), sep = "_") %>%
          mutate(set = paste0(set, " ", set_value)) %>%
          dplyr::select(-c(matrix, row, column, set_value)) %>%
          separate(col = scheme_value, into = c("sol_type", "dil", "rep_bio", "name"), sep = ",")

        t5 <- extract_matrix_values("Time 5_Set", input$n_col, input$n_row, input) %>%
          mutate(id = paste0(row, column)) %>%
          rename(biolu = value) %>%
          dplyr::select(-c(row, column)) %>%
          separate(col = matrix, into = c("time", "set", "set_value"), sep = "_") %>%
          mutate(set = paste0(set, " ", set_value)) %>%
          dplyr::select(-c(set_value))

        t15 <- extract_matrix_values("Time 15_Set", input$n_col, input$n_row, input) %>%
          mutate(id = paste0(row, column)) %>%
          rename(biolu = value) %>%
          dplyr::select(-c(row, column)) %>%
          separate(col = matrix, into = c("time", "set", "set_value"), sep = "_") %>%
          mutate(set = paste0(set, " ", set_value)) %>%
          dplyr::select(-c(set_value))

        t30 <- extract_matrix_values("Time 30_Set", input$n_col, input$n_row, input) %>%
          mutate(id = paste0(row, column)) %>%
          rename(biolu = value) %>%
          dplyr::select(-c(row, column)) %>%
          separate(col = matrix, into = c("time", "set", "set_value"), sep = "_") %>%
          mutate(set = paste0(set, " ", set_value)) %>%
          dplyr::select(-c(set_value))

        data_5 <- merge(scheme, t5, by = c("id", "set"), all = TRUE)
        data_15 <- merge(scheme, t15, by = c("id", "set"), all = TRUE)
        data_30 <- merge(scheme, t30, by = c("id", "set"), all = TRUE)

        final_data <- rbind(data_5, data_15, data_30)
        final_data$manip_name <- input$manip_name
        final_data_corrected <- correct_final_data(final_data, input$neg_control)
        final_data_corrected <- final_data_corrected %>% filter(!is.na(biolu))
        final_data_corrected <- final_data_corrected %>% dplyr::select(c(manip_name, set, time, id, sol_type, name, dil, rep_bio,
                                                                         biolu, mean_biolu_nacl, biolu_corr, perc_inhib, perc_inhib_corr))
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
      content = function(file) {write.csv(final_data_corrected(), file, row.names = FALSE)}
    )



    ### e. Création d'un graphique final "général" ####



    # Adapter les données pour le graphique final
    final_data_corrected <-

      reactive(

        {

          extract_matrix_values <-

            function(matrix_prefix, ncol, nrow, input) {

              matrix_data <- data.frame(
                matrix = character(),
                row = character(),
                column = integer(),
                value = character(),
                stringsAsFactors = FALSE
              )

              for (i in 1:input$n_set){
                matrix_name <- paste0(matrix_prefix, "_", i)
                for (row in LETTERS[1:nrow]){
                  for (col in 1:ncol){
                    cell_id <- paste(matrix_name, row, col, sep = "_")
                    value <- input[[cell_id]]
                    matrix_data <- rbind(matrix_data, data.frame(
                      matrix = matrix_name,
                      row = row,
                      column = col,
                      value = value,
                      stringsAsFactors = FALSE)
                    )

                  }
                }
              }

              return(matrix_data)

            }

          observe({
            sol_types <- unique(final_data_corrected()$sol_type)
            updateSelectInput(session, "sol_type", choices = sol_types)
          })

          scheme <- extract_matrix_values("Scheme_Set", input$n_col, input$n_row, input) %>%
            mutate(id = paste0(row, column)) %>%
            rename(scheme_value = value) %>%
            separate(col = matrix, into = c("matrix", "set", "set_value"), sep = "_") %>%
            mutate(set = paste0(set, " ", set_value)) %>%
            dplyr::select(-c(matrix, row, column, set_value)) %>%
            separate(col = scheme_value, into = c("sol_type", "dil", "rep_bio", "name"), sep = ",")

          t5 <- extract_matrix_values("Time 5_Set", input$n_col, input$n_row, input) %>%
            mutate(id = paste0(row, column)) %>%
            rename(biolu = value) %>%
            dplyr::select(-c(row, column)) %>%
            separate(col = matrix, into = c("time", "set", "set_value"), sep = "_") %>%
            mutate(set = paste0(set, " ", set_value)) %>%
            dplyr::select(-c(set_value))

          t15 <- extract_matrix_values("Time 15_Set", input$n_col, input$n_row, input) %>%
            mutate(id = paste0(row, column)) %>%
            rename(biolu = value) %>%
            dplyr::select(-c(row, column)) %>%
            separate(col = matrix, into = c("time", "set", "set_value"), sep = "_") %>%
            mutate(set = paste0(set, " ", set_value)) %>%
            dplyr::select(-c(set_value))

          t30 <- extract_matrix_values("Time 30_Set", input$n_col, input$n_row, input) %>%
            mutate(id = paste0(row, column)) %>%
            rename(biolu = value) %>%
            dplyr::select(-c(row, column)) %>%
            separate(col = matrix, into = c("time", "set", "set_value"), sep = "_") %>%
            mutate(set = paste0(set, " ", set_value)) %>%
            dplyr::select(-c(set_value))

          data_5 <- merge(scheme, t5, by = c("id", "set"), all = TRUE)
          data_15 <- merge(scheme, t15, by = c("id", "set"), all = TRUE)
          data_30 <- merge(scheme, t30, by = c("id", "set"), all = TRUE)

          final_data <- rbind(data_5, data_15, data_30)
          final_data$manip_name <- input$manip_name
          final_data_corrected <- correct_final_data(final_data, input$neg_control)
          final_data_corrected <- final_data_corrected %>% filter(!is.na(biolu))
          final_data_corrected <- final_data_corrected %>% dplyr::select(c(manip_name, set, time, id, sol_type, name, dil, rep_bio,
                                                                           biolu, mean_biolu_nacl, biolu_corr, perc_inhib, perc_inhib_corr))
          final_data_corrected <- final_data_corrected %>% arrange(dil)

          niveaux <- sort(as.numeric(unique(final_data_corrected$dil)))

          return(final_data_corrected)

        })



    # Sélectionner certaines variables par défaut pour le graphique
    observe({
      updateSelectInput(session, "x_variable", choices = c("sol_type", "dil", "rep_bio", "name", "id", "time", "set", "manip_name"), selected = "sol_type")
      updateSelectInput(session, "y_variable", choices = c("biolu", "mean_biolu_nacl", "biolu_corr", "perc_inhib", "perc_inhib_corr"), selected = "perc_inhib_corr")
      updateSelectInput(session, "color_variable", choices = c("sol_type", "dil", "rep_bio", "name", "id", "time", "set", "manip_name"), selected = "dil")
    })



    # Choix des substances dynamique en fonction des valeurs entrées par l'utilisateur
    output$dynamic_substance_checkboxes <-

      renderUI({

        unique_substances <- unique(final_data_corrected()$sol_type)

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
          label = "Dilutions to include:",
          choices = unique_dilutions,
          selected = unique_dilutions
        )

        return(checkbox_group)

      })



    # Adapter les données en fonction des choix pour le graphique
    correct_data_microtox_filtered <-

      reactive({

        final_data_corrected <- final_data_corrected()
        final_data_corrected <- filter_data(final_data_corrected, input$filter_t5, input$filter_t15, input$filter_t30)
        selected_dilutions <- as.numeric(input$selected_dilutions)
        final_data_corrected <- filter_data_dil(final_data_corrected, selected_dilutions)

        return(final_data_corrected)

      })



    # Création du graphique final
    output$plot_biolutoxR <-

      renderPlot({

        req(input$x_variable, input$y_variable, input$color_variable)

        data_plot <- correct_data_microtox_filtered()
        data_plot$dil = factor(data_plot$dil, levels = sort(as.numeric(unique(data_plot$dil))))
        data_plot$time = factor(data_plot$time, levels = c("Time 5", "Time 15", "Time 30"))

        selected_substances <- input$selected_substances

        if (!is.null(selected_substances) && length(selected_substances) > 0) {
          data_plot <- data_plot %>% filter(sol_type %in% selected_substances)
        }

        p_data <- ggplot(data_plot, aes_string(x = input$x_variable,
                                               y = input$y_variable,
                                               color = input$color_variable)) +
          geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.85)) +
          geom_jitter(alpha = 0.5, position = position_dodge(width = 0.85)) +
          stat_summary(fun = "mean", geom = "point", shape = 4, size = 3, position = position_dodge(width = 0.85)) +
          theme_test() +
          theme(text = element_text(size = 20))
        print(p_data)

      })



    ### f. Création d'un graphique final "toxicité" ####



    # Création du graphique final
    output$ec50_plot <-

      renderPlot({

        tryCatch({
          ec_data <- final_data_corrected()
          ec_data <- filter_data(ec_data, input$filter_t5_2, input$filter_t15_2, input$filter_t30_2)
          table_ec <- toxicity_data_table(ec_data, substance = input$sol_type)
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
        ec_data <- filter_data(ec_data, input$filter_t5_2, input$filter_t15_2, input$filter_t30_2)
        table_ec <- toxicity_data_table(ec_data, substance = input$sol_type)
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
        ec_data <- filter_data(ec_data, input$filter_t5_2, input$filter_t15_2, input$filter_t30_2)
        table_ec <- toxicity_data_table(ec_data, substance = input$sol_type)
        ec_X_value <- toxicity_data_ecX(table_ec, X = input$ecX_input_value)
        paste0("More specifically, bacterial bioluminescence was inhibited for ",  input$ecX_input_value, "% of organisms the batch exposed to the ", input$sol_type, " at a dilution of ", ec_X_value, ".")
      }, error = function(e) {
        "Changer de substance"
      })

    })



    ### g. Création d'un graphique final "toxicité" ####
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
# ------------------------ LANCEMENT example.biolutoxR() -----------------------
# ############################################################################ #



# Run la function example.biolutoxR()
example.biolutoxR()


