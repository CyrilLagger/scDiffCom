ui <- function(
  request
) {
  tagList(
    fluidPage(
      theme = shinythemes::shinytheme("cerulean"),
      titlePanel(
        title = tags$table(
          style = "width: 100%",
          tags$tbody(
            tags$tr(
              tags$td(
                tags$span(
                  style = "font-size: 26px",
                  paste(
                    "Shiny Report for scDiffCom object",
                    .object_@parameters$object_name
                  )
                )
              )
            )
          )
        ),
        windowTitle = "ShinyReport"
      ),
      navbarPage(
        title = NULL,
        id = "navbarID",
        tabPanel(
          title = "Results",
          mod_results_ui("results_ui_1"),
          value = "results_navbar"
        ),
        tabPanel(
          title = "Parameters",
          mod_parameters_ui("parameters_ui_1"),
          value = "parameters_navbar"
        )
      )
    )
  )
}

mod_parameters_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(
        column(
          width = 8,
          offset = 2,
          DT::DTOutput(
            outputId = ns("PARAMETERS_DT")
          )
        )
      )
    )
  )
}

mod_results_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(
        column(
          style = "text-align:center;",
          width = 6,
          titlePanel(htmlOutput(ns("RESULTS_TITLE"))),
          offset = 3
        ),
      )
    ),
    tabsetPanel(
      type = "tabs",
      mod_results_cci_ui(ns("results_cci_ui_1")),
      mod_results_ora_ui(ns("results_ora_ui_1")),
      id = "active_results_panel",
      selected = "RESULTS_INTERACTION_ANALYSIS"
    )
  )
}

mod_results_cci_ui <- function(id){
  ns <- NS(id)
  tabPanel(
    title = "Detected Interactions",
    sidebarLayout(
      sidebarPanel(
        width = 2,
        downloadButton(
          ns("RESULTS_DOWNLOAD_TABLE"),
          "Download Full Table"
        ),
        hr(
          style = "border-top: 1px solid #000000;"
        ),
        h4("Filtering Options"),
        uiOutput(ns("RESULTS_EMITTER_CHOICE")),
        uiOutput(ns("RESULTS_RECEIVER_CHOICE")),
        selectizeInput(
          ns("RESULTS_LRI_CHOICE"),
          choices = NULL,
          label = "Ligand-Receptor Interactions",
          multiple = TRUE
        ),
        selectizeInput(
          ns("RESULTS_GENE_CHOICE"),
          choices = NULL,
          label = "Individual Genes",
          multiple = TRUE
        ),
        selectizeInput(
          ns("RESULTS_GO_CHOICE"),
          choices = NULL,
          label = "GO Terms",
          multiple = TRUE
        ),
        selectizeInput(
          ns("RESULTS_KEGG_CHOICE"),
          choices = NULL,
          label = "KEGG Pathways",
          multiple = TRUE
        ),
        actionButton(
          inputId = ns("RESULTS_FILTER_BUTTON"),
          label = "Filter"
        ),
        actionButton(
          inputId = ns("RESULTS_RESET_BUTTON"),
          label = "Undo Filtering"
        )
      ),
      mainPanel(
        width = 10,
        uiOutput(ns("RESULTS_CCI_TITLE")),
        uiOutput(ns("RESULTS_CCI_DETAILS"))
      )
    ),
    value = "RESULTS_INTERACTION_ANALYSIS"
  )
}

mod_results_ora_ui <- function(id){
  ns <- NS(id)
  tabPanel(
    title = "Over-representation Analysis",
    sidebarLayout(
      sidebarPanel(
        width = 2,
        selectInput(
          inputId = ns("RESULTS_ORA_CATEGORY_CHOICE"),
          label = "Category",
          choices = c(
            "By Cell Types",
            "By GO/KEGG",
            "By Genes"
          )
        ),
        hr(
          style = "border-top: 1px solid #000000;"
        ),
        selectInput(
          inputId = ns("RESULTS_ORA_TYPE_CHOICE"),
          label = "Age Regulation",
          choices = c("UP", "DOWN", "FLAT")
        )
      ),
      mainPanel(
        width = 10,
        uiOutput(ns("RESULTS_ORA_TITLE")),
        uiOutput(ns("RESULTS_ORA_DETAILS"))
      )
    ),
    value = "RESULTS_ORA"
  )
}
