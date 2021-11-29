server <- function(
  input,
  output,
  session
) {
  mod_results_server("results_ui_1")
  mod_parameters_server("parameters_ui_1")
}

mod_parameters_server <- function(id) {
  moduleServer(
    id,
    function(
      input,
      output,
      session
    ) {
      ns <- session$ns

      output$PARAMETERS_DT <- DT::renderDT({
        display_parameters_table(.object_@parameters)
      })

    })
}

mod_results_server <- function(id) {
  moduleServer(
    id,
    function(
      input,
      output,
      session
    ) {
      ns <- session$ns

      output$RESULTS_TITLE <- renderUI({
        tags$div(
          style = paste(
            #"display: inline-block;",
            #"text-align: center;",
            "font-size: 26px"
          ),
          "Results"
        )
      })

      mod_results_cci_server("results_cci_ui_1")
      mod_results_ora_server("results_ora_ui_1")

    })
}

mod_results_cci_server <- function(id, rv_hidden) {
  moduleServer(
    id,
    function(
      input,
      output,
      session
    ) {
      ns <- session$ns

      output$RESULTS_DOWNLOAD_TABLE <- downloadHandler(
        filename = function() {
          paste0(
            "cci_table_",
            tolower(
              gsub(
                " ",
                "_",
                .object_@parameters$object_name,
                fixed = TRUE
              )
            ),
            ".csv"
          )
        },
        content = function(file) {
          CCI_table_downl <- GetTableCCI(
            object = .object_,
            type = "detected",
            simplified = TRUE
          )
          fwrite(CCI_table_downl, file)
        }
      )

      output$RESULTS_EMITTER_CHOICE <- renderUI({
        choices <- sort(
          unique(
            .object_@cci_table_detected$EMITTER_CELLTYPE
          )
        )
        shinyWidgets::pickerInput(
          inputId = ns("RESULTS_EMITTER_CHOICE"),
          label = "Emitter Cell Types",
          choices = choices,
          selected = choices,
          options = list(`actions-box` = TRUE),
          multiple = TRUE
        )
      })

      output$RESULTS_RECEIVER_CHOICE <- renderUI({
        choices <- sort(
          unique(
            .object_@cci_table_detected$RECEIVER_CELLTYPE
          )
        )
        shinyWidgets::pickerInput(
          inputId = ns("RESULTS_RECEIVER_CHOICE"),
          label = "Receiver Cell Types",
          choices = choices,
          selected = choices,
          options = list(`actions-box` = TRUE),
          multiple = TRUE
        )
      })

      output$RESULTS_LRI_CHOICE <- renderUI({
        ALL_LRI_LABEL <- "All LRIs"
        choices <-
          c(
            ALL_LRI_LABEL,
            sort(unique(.object_@cci_table_detected$LRI))
          )
        selectizeInput(
          inputId = ns("RESULTS_LRI_CHOICE"),
          label = "Ligand-Receptor Interactions",
          choices = choices,
          selected = ALL_LRI_LABEL,
          multiple = TRUE,
          options = list(
            allowEmptyOption = TRUE,
            placeholder = 'Type LRIs',
            maxOptions = length(choices)
          )
        )
      })

      output$RESULTS_GENE_CHOICE <- renderUI({
        ALL_GENE_LABEL <- "All Genes"
        choices <- c(
          ALL_GENE_LABEL,
          sort(
            unique(
              c(
                sapply(
                  seq_along(.object_@parameters$max_nL),
                  function(i) {
                    .object_@cci_table_detected[[paste0("LIGAND_", i)]]
                  }
                ),
                sapply(
                  seq_along(.object_@parameters$max_nR),
                  function(i) {
                    .object_@cci_table_detected[[paste0("RECEPTOR_", i)]]
                  }
                )
              )
            )
          )
        )
        selectizeInput(
          inputId = ns("RESULTS_GENE_CHOICE"),
          label = "Individual Genes",
          choices = choices,
          selected = ALL_GENE_LABEL,
          multiple = TRUE,
          options = list(
            allowEmptyOption = TRUE,
            placeholder = 'Type Genes',
            maxOptions = length(choices)
          )
        )
      })

      output$RESULTS_GO_CHOICE <- renderUI({
        ALL_GO_LABEL <- "All GO Terms"
        choices_go <- c(
          ALL_GO_LABEL,
          sort(unique(.object_@ora_table$GO_TERMS$VALUE))
        )
        selectizeInput(
          inputId = ns("RESULTS_GO_CHOICE"),
          label = "GO Terms",
          choices = choices_go,
          selected = ALL_GO_LABEL,
          multiple = TRUE,
          options = list(
            allowEmptyOption = TRUE,
            placeholder = "Type GO Terms",
            maxOptions = length(choices_go)
          )
        )
      })

      output$RESULTS_KEGG_CHOICE <- renderUI({
        ALL_KEGG_LABEL = 'All KEGG Pathways'
        choices <- c(
          ALL_KEGG_LABEL,
          sort(unique(.object_@ora_table$KEGG_PWS$VALUE))
        )
        selectizeInput(
          inputId = ns("RESULTS_KEGG_CHOICE"),
          label = "KEGG Pathways",
          choices = choices,
          selected = ALL_KEGG_LABEL,
          multiple = TRUE,
          options = list(
            allowEmptyOption = TRUE,
            placeholder = 'Type KEGG Pathways',
            maxOptions = length(choices)
          )
        )
      })

      output$RESULTS_CCI_TITLE <- renderUI({
        fluidPage(
          fluidRow(
            column(
              width = 12,
              titlePanel(
                tags$p(
                  div(
                    style = paste(
                      "width: 80%;",
                      "margin:auto;",
                      "font-size: 20px;",
                      "text-align: center;",
                      "margin-bottom: 50px;"
                    ),
                    "Plots and Table"
                  )
                )
              )
            )
          )
        )
      })

      output$RESULTS_CCI_DETAILS <- renderUI({
        fluidPage(
          fluidRow(
            column(
              style = "padding: 10px;margin-bottom:50px;",
              width = 5,
              offset = 1,
              plotly::plotlyOutput(
                outputId = ns("RESULTS_PLOTLY_VOLCANO"),
                height = "460px"
              )
            ),
            column(
              style = "padding: 10px;margin-bottom:50px;",
              width = 5,
              offset = 1,
              plotly::plotlyOutput(
                outputId = ns("RESULTS_PLOTLY_SCORE"),
                height = "460px"
              )
            )
          ),
          fluidRow(
            column(
              style = "padding: 10px;",
              width = 5,
              offset = 3,
              plotly::plotlyOutput(
                outputId = ns("RESULTS_PLOTLY_LRFC"),
                height = "460px"
              )
            )
          ),
          fluidRow(
            column(
              width = 10,
              offset = 1,
              DT::DTOutput(
                outputId = ns("RESULTS_CCI_DT")
              )
            )
          )
        )
      })

      output$RESULTS_PLOTLY_VOLCANO <- plotly::renderPlotly({
        plot_volcano_CCI(CCI_table(), .object_@parameters)
      })

      output$RESULTS_PLOTLY_SCORE <- plotly::renderPlotly({
        plot_scores_CCI(CCI_table(), .object_@parameters)
      })

      output$RESULTS_PLOTLY_LRFC <- plotly::renderPlotly({
        plot_lrfc_CCI(CCI_table(), .object_@parameters)
      })

      output$RESULTS_CCI_DT <- DT::renderDT({
        display_CCI_table(CCI_table(), .object_@parameters)
      })

      CCI_table <- reactive({
        if (filter_values$do_filtering) {
          CCI_table <- subset_CCI_table(
            CCI_table = .object_@cci_table_detected,
            emitter_choice = filter_values$emitter_choice,
            receiver_choice = filter_values$receiver_choice,
            LRI_choice = filter_values$LRI_choice,
            GENE_choice = filter_values$GENE_choice,
            GO_choice = filter_values$GO_choice,
            KEGG_choice = filter_values$KEGG_choice,
            filter = TRUE
          )
        } else {
          CCI_table <- subset_CCI_table(
            CCI_table = .object_@cci_table_detected,
            filter = FALSE
          )
        }
        CCI_table
      })

      filter_values <- reactiveValues(
        do_filtering = FALSE,
        emitter_choice = NULL,
        receiver_choice = NULL,
        LRI_choice = NULL,
        GENE_choice = NULL,
        GO_choice = NULL,
        KEGG_choice = NULL
      )

      observeEvent(
        input$RESULTS_FILTER_BUTTON,
        {
          filter_values$do_filtering <- TRUE
          filter_values$emitter_choice <- input$RESULTS_EMITTER_CHOICE
          filter_values$receiver_choice <- input$RESULTS_RECEIVER_CHOICE
          filter_values$LRI_choice <- input$RESULTS_LRI_CHOICE
          filter_values$GENE_choice <- input$RESULTS_GENE_CHOICE
          filter_values$GO_choice <- input$RESULTS_GO_CHOICE
          filter_values$KEGG_choice <- input$RESULTS_KEGG_CHOICE
        }
      )

      observeEvent(
        input$RESULTS_RESET_BUTTON,
        {
          filter_values$do_filtering <- FALSE
          filter_values$emitter_choice <- NULL
          filter_values$receiver_choice <- NULL
          filter_values$LRI_choice <- NULL
          filter_values$GENE_choice <- NULL
          filter_values$GO_choice <- NULL
          filter_values$KEGG_choice <- NULL
          choices_cts <- sort(
            unique(
              .object_@cci_table_detected$EMITTER_CELLTYPE
            )
          )
          shinyWidgets::updatePickerInput(
            session = session,
            inputId = "RESULTS_EMITTER_CHOICE",
            choices = choices_cts,
            selected = choices_cts
          )
          shinyWidgets::updatePickerInput(
            session = session,
            inputId = "RESULTS_RECEIVER_CHOICE",
            choices = choices_cts,
            selected = choices_cts
          )
          ALL_LRI_LABEL <- "All LRIs"
          choices_lri <-
            c(
              ALL_LRI_LABEL,
              sort(unique(.object_@cci_table_detected$LRI))
            )
          updateSelectizeInput(
            session = session,
            "RESULTS_LRI_CHOICE",
            choices = choices_lri,
            selected = ALL_LRI_LABEL,
            options = list(
              allowEmptyOption = TRUE,
              placeholder = "Type LRIs",
              maxOptions = length(choices_lri)
            ),
            server = TRUE
          )
          ALL_GENE_LABEL <- "All Genes"
          choices_gene <- c(
            ALL_GENE_LABEL,
            sort(
              unique(
                c(
                  sapply(
                    seq_along(.object_@parameters$max_nL),
                    function(i) {
                      .object_@cci_table_detected[[paste0("LIGAND_", i)]]
                    }
                  ),
                  sapply(
                    seq_along(.object_@parameters$max_nR),
                    function(i) {
                      .object_@cci_table_detected[[paste0("RECEPTOR_", i)]]
                    }
                  )
                )
              )
            )
          )
          updateSelectizeInput(
            session = session,
            "RESULTS_GENE_CHOICE",
            choices = choices_gene,
            selected = ALL_GENE_LABEL,
            options = list(
              allowEmptyOption = TRUE,
              placeholder = "Type Genes",
              maxOptions = length(choices_gene)
            ),
            server = TRUE
          )
          ALL_GO_LABEL <- "All GO Terms"
          choices_go <- c(
            ALL_GO_LABEL,
            sort(unique(.object_@ora_table$GO_TERMS$VALUE))
          )
          updateSelectizeInput(
            session = session,
            "RESULTS_GO_CHOICE",
            choices = choices_go,
            selected = ALL_GO_LABEL,
            options = list(
              allowEmptyOption = TRUE,
              placeholder = "Type GO Terms",
              maxOptions = length(choices_go)
            ),
            server = TRUE
          )
          ALL_KEGG_LABEL = "All KEGG Pathways"
          choices_kegg <- c(
            ALL_KEGG_LABEL,
            sort(unique(.object_@ora_table$KEGG_PWS$VALUE))
          )
          updateSelectizeInput(
            session = session,
            "RESULTS_KEGG_CHOICE",
            choices = choices_kegg,
            selected = ALL_KEGG_LABEL,
            options = list(
              allowEmptyOption = TRUE,
              placeholder = "Type KEGG Pathways",
              maxOptions = length(choices_kegg)
            ),
            server = TRUE
          )
        }
      )
    })
}

mod_results_ora_server <- function(id) {
  moduleServer(
    id,
    function(
      input,
      output,
      session
    ) {
      ns <- session$ns

      output$RESULTS_ORA_TITLE <- renderUI({
        fluidPage(
          fluidRow(
            column(
              width = 12,
              titlePanel(
                tags$p(
                  div(
                    style = paste(
                      "width: 80%;",
                      "margin:auto;",
                      "font-size: 20px;",
                      "text-align: center;",
                      "margin-bottom:50px"
                    ),
                    "Over-representation Results",
                  )
                )
              )
            )
          )
        )
      })

      output$RESULTS_ORA_DETAILS <-  renderUI({
        if (input$RESULTS_ORA_CATEGORY_CHOICE == "By Cell Types") {
          fluidPage(
            fluidRow(
              column(
                width = 9,
                offset = 1,
                style = "margin-bottom: 50px;",
                visNetwork::visNetworkOutput(
                  ns("RESULTS_ORA_NETWORK_PLOT"),
                  height = "700px"
                )
              )
            ),
            fluidRow(
              column(
                style = "padding: 10px;",
                width = 8,
                offset = 2,
                DT::DTOutput(ns("RESULTS_ORA_TABLE_ERI"))
              )
            ),
            fluidRow(
              column(
                style = "padding: 10px;",
                width = 8,
                offset = 2,
                DT::DTOutput(ns("RESULTS_ORA_TABLE_EMITTER"))
              )
            ),
            fluidRow(
              column(
                style = "padding: 10px;",
                width = 8,
                offset = 2,
                DT::DTOutput(ns("RESULTS_ORA_TABLE_RECEIVER"))
              )
            )
          )
        } else if (input$RESULTS_ORA_CATEGORY_CHOICE == "By Genes") {
          fluidPage(
            fluidRow(
              style = "margin-bottom:50px;",
              column(
                width = 4,
                plotOutput(
                  ns("RESULTS_ORA_PLOT_LRI"),
                  height = "500px"
                )
              ),
              column(
                width = 4,
                plotOutput(
                  ns("RESULTS_ORA_PLOT_LIGAND"),
                  height = "500px"
                )
              ),
              column(
                width = 4,
                plotOutput(
                  ns("RESULTS_ORA_PLOT_RECEPTOR"),
                  height = "500px"
                )
              )
            ),
            fluidRow(
              column(
                style = "padding: 10px;",
                width = 8,
                offset = 2,
                DT::DTOutput(ns("RESULTS_ORA_TABLE_LRI"))
              )
            ),
            fluidRow(
              column(
                style = "padding: 10px;",
                width = 8,
                offset = 2,
                DT::DTOutput(ns("RESULTS_ORA_TABLE_LIGAND"))
              )
            ),
            fluidRow(
              column(
                style = "padding: 10px;",
                width = 8,
                offset = 2,
                DT::DTOutput(ns("RESULTS_ORA_TABLE_RECEPTOR"))
              )
            )
          )
        } else if (input$RESULTS_ORA_CATEGORY_CHOICE == "By GO/KEGG") {
          if (.reduce_go_) {
            fluidPage(
              fluidRow(
                style = "margin-bottom:50px;",
                column(
                  width = 6,
                  plotly::plotlyOutput(
                    ns("RESULTS_ORA_TREEMAP_GO_BP"),
                    height = "700px"
                  )
                ),
                column(
                  width = 6,
                  plotly::plotlyOutput(
                    ns("RESULTS_ORA_TREEMAP_GO_MF"),
                    height = "700px"
                  )
                )
              ),
              fluidRow(
                column(
                  width = 6,
                  plotly::plotlyOutput(
                    ns("RESULTS_ORA_TREEMAP_GO_CC"),
                    height = "700px"
                  )
                ),
                column(
                  width = 6,
                  plotOutput(
                    ns("RESULTS_ORA_PLOT_KEGG"),
                    height = "700px")
                )
              ),
              fluidRow(
                column(
                  style = "padding: 10px;",
                  width = 8,
                  offset = 2,
                  DT::DTOutput(
                    ns("RESULTS_ORA_TABLE_GO")
                  )
                )
              ),
              fluidRow(
                column(
                  style = "padding: 10px;",
                  width = 8,
                  offset = 2,
                  DT::DTOutput(ns("RESULTS_ORA_TABLE_KEGG"))
                )
              )
            )
          } else {
            fluidPage(
              fluidRow(
                style = "margin-bottom:50px;",
                column(
                  width = 6,
                  plotOutput(
                    ns("RESULTS_ORA_PLOT_GO_BP"),
                    height = "700px"
                  )
                ),
                column(
                  width = 6,
                  plotOutput(
                    ns("RESULTS_ORA_PLOT_GO_MF"),
                    height = "700px"
                  )
                )
              ),
              fluidRow(
                column(
                  width = 6,
                  plotOutput(
                    ns("RESULTS_ORA_PLOT_GO_CC"),
                    height = "700px"
                  )
                ),
                column(
                  width = 6,
                  plotOutput(
                    ns("RESULTS_ORA_PLOT_KEGG"),
                    height = "700px")
                )
              ),
              fluidRow(
                column(
                  style = "padding: 10px;",
                  width = 8,
                  offset = 2,
                  DT::DTOutput(
                    ns("RESULTS_ORA_TABLE_GO")
                  )
                )
              ),
              fluidRow(
                column(
                  style = "padding: 10px;",
                  width = 8,
                  offset = 2,
                  DT::DTOutput(ns("RESULTS_ORA_TABLE_KEGG"))
                )
              )
            )
          }
        }
      })

      output$RESULTS_ORA_NETWORK_PLOT <- visNetwork::renderVisNetwork({
        scDiffCom::BuildNetwork(
          object = .object_
        )
      })

      output$RESULTS_ORA_TABLE_ERI <- DT::renderDT({
        req(
          input$RESULTS_ORA_CATEGORY_CHOICE,
          input$RESULTS_ORA_TYPE_CHOICE
        )
        display_ORA_table(
          ORA_table = .object_@ora_table$ER_CELLTYPES,
          category_choice = "ER_CELLTYPES",
          type_choice = input$RESULTS_ORA_TYPE_CHOICE
        )
      })

      output$RESULTS_ORA_TABLE_EMITTER <- DT::renderDT({
        req(
          input$RESULTS_ORA_CATEGORY_CHOICE,
          input$RESULTS_ORA_TYPE_CHOICE
        )
        display_ORA_table(
          ORA_table = .object_@ora_table$EMITTER_CELLTYPE,
          category_choice = "EMITTER_CELLTYPE",
          type_choice = input$RESULTS_ORA_TYPE_CHOICE
        )
      })

      output$RESULTS_ORA_TABLE_RECEIVER <- DT::renderDT({
        req(
          input$RESULTS_ORA_CATEGORY_CHOICE,
          input$RESULTS_ORA_TYPE_CHOICE
        )
        display_ORA_table(
          ORA_table = .object_@ora_table$RECEIVER_CELLTYPE,
          category_choice = "RECEIVER_CELLTYPE",
          type_choice = input$RESULTS_ORA_TYPE_CHOICE
        )
      })

      output$RESULTS_ORA_PLOT_LRI <- renderPlot({
        req(
          input$RESULTS_ORA_CATEGORY_CHOICE,
          input$RESULTS_ORA_TYPE_CHOICE
        )
        scDiffCom::PlotORA(
          object = .object_,
          category = "LRI",
          regulation = input$RESULTS_ORA_TYPE_CHOICE,max_terms_show = 20,
          OR_threshold = 1,
          bh_p_value_threshold = 0.05
        )
      })

      output$RESULTS_ORA_PLOT_LIGAND <- renderPlot({
        req(
          input$RESULTS_ORA_CATEGORY_CHOICE,
          input$RESULTS_ORA_TYPE_CHOICE
        )
        scDiffCom::PlotORA(
          object = .object_,
          category = "LIGAND_COMPLEX",
          regulation = input$RESULTS_ORA_TYPE_CHOICE,max_terms_show = 20,
          OR_threshold = 1,
          bh_p_value_threshold = 0.05
        )
      })

      output$RESULTS_ORA_PLOT_RECEPTOR <- renderPlot({
        req(
          input$RESULTS_ORA_CATEGORY_CHOICE,
          input$RESULTS_ORA_TYPE_CHOICE
        )
        scDiffCom::PlotORA(
          object = .object_,
          category = "RECEPTOR_COMPLEX",
          regulation = input$RESULTS_ORA_TYPE_CHOICE,max_terms_show = 20,
          OR_threshold = 1,
          bh_p_value_threshold = 0.05
        )
      })

      output$RESULTS_ORA_TABLE_LRI <- DT::renderDT({
        req(
          input$RESULTS_ORA_CATEGORY_CHOICE,
          input$RESULTS_ORA_TYPE_CHOICE
        )
        display_ORA_table(
          ORA_table = .object_@ora_table$LRI,
          category_choice = "LRI",
          type_choice = input$RESULTS_ORA_TYPE_CHOICE
        )
      })

      output$RESULTS_ORA_TABLE_LIGAND <- DT::renderDT({
        req(
          input$RESULTS_ORA_CATEGORY_CHOICE,
          input$RESULTS_ORA_TYPE_CHOICE
        )
        display_ORA_table(
          ORA_table = .object_@ora_table$LIGAND_COMPLEX,
          category_choice = "LIGAND_COMPLEX",
          type_choice = input$RESULTS_ORA_TYPE_CHOICE
        )
      })

      output$RESULTS_ORA_TABLE_RECEPTOR <- DT::renderDT({
        req(
          input$RESULTS_ORA_CATEGORY_CHOICE,
          input$RESULTS_ORA_TYPE_CHOICE
        )
        display_ORA_table(
          ORA_table = .object_@ora_table$RECEPTOR_COMPLEX,
          category_choice = "RECEPTOR_COMPLEX",
          type_choice = input$RESULTS_ORA_TYPE_CHOICE
        )
      })

      output$RESULTS_ORA_PLOT_GO_BP <- renderPlot({
        req(
          input$RESULTS_ORA_CATEGORY_CHOICE,
          input$RESULTS_ORA_TYPE_CHOICE
        )
        scDiffCom::PlotORA(
          object = .object_,
          category = "GO_TERMS",
          regulation = input$RESULTS_ORA_TYPE_CHOICE,
          max_terms_show = 20,
          GO_aspect = "biological_process",
          OR_threshold = 1,
          bh_p_value_threshold = 0.05
        )
      })

      output$RESULTS_ORA_PLOT_GO_MF <- renderPlot({
        req(
          input$RESULTS_ORA_CATEGORY_CHOICE,
          input$RESULTS_ORA_TYPE_CHOICE
        )
        scDiffCom::PlotORA(
          object = .object_,
          category = "GO_TERMS",
          regulation = input$RESULTS_ORA_TYPE_CHOICE,
          max_terms_show = 20,
          GO_aspect = "molecular_function",
          OR_threshold = 1,
          bh_p_value_threshold = 0.05
        )
      })

      output$RESULTS_ORA_PLOT_GO_CC <- renderPlot({
        req(
          input$RESULTS_ORA_CATEGORY_CHOICE,
          input$RESULTS_ORA_TYPE_CHOICE
        )
        scDiffCom::PlotORA(
          object = .object_,
          category = "GO_TERMS",
          regulation = input$RESULTS_ORA_TYPE_CHOICE,
          max_terms_show = 20,
          GO_aspect = "cellular_component",
          OR_threshold = 1,
          bh_p_value_threshold = 0.05
        )
      })

      output$RESULTS_ORA_TREEMAP_GO_BP <- plotly::renderPlotly({
        req(
          input$RESULTS_ORA_CATEGORY_CHOICE,
          input$RESULTS_ORA_TYPE_CHOICE
        )
        plot_ORA_GO_treemap(
          GO_REDUCED_table = .reduced_go_table_,
          type_choice = input$RESULTS_ORA_TYPE_CHOICE,
          go_aspect_choice = "biological_process",
          title_text = paste0(
            "GO Biological Processes - ",
            input$RESULTS_ORA_TYPE_CHOICE
          )
        )
      })

      output$RESULTS_ORA_TREEMAP_GO_MF <- plotly::renderPlotly({
        req(
          input$RESULTS_ORA_CATEGORY_CHOICE,
          input$RESULTS_ORA_TYPE_CHOICE
        )
        plot_ORA_GO_treemap(
          GO_REDUCED_table = .reduced_go_table_,
          type_choice = input$RESULTS_ORA_TYPE_CHOICE,
          go_aspect_choice = "molecular_function",
          title_text = paste0(
            "GO Molecular Functions - ",
            input$RESULTS_ORA_TYPE_CHOICE
          )
        )
      })

      output$RESULTS_ORA_TREEMAP_GO_CC <- plotly::renderPlotly({
        req(
          input$RESULTS_ORA_CATEGORY_CHOICE,
          input$RESULTS_ORA_TYPE_CHOICE
        )
        plot_ORA_GO_treemap(
          GO_REDUCED_table = .reduced_go_table_,
          type_choice = input$RESULTS_ORA_TYPE_CHOICE,
          go_aspect_choice = "cellular_component",
          title_text = paste0(
            "GO Cellular Components - ",
            input$RESULTS_ORA_TYPE_CHOICE
          )
        )
      })

      output$RESULTS_ORA_PLOT_KEGG <- renderPlot({
        req(
          input$RESULTS_ORA_CATEGORY_CHOICE,
          input$RESULTS_ORA_TYPE_CHOICE
        )
        scDiffCom::PlotORA(
          object = .object_,
          category = "KEGG_PWS",
          regulation = input$RESULTS_ORA_TYPE_CHOICE,max_terms_show = 20,
          OR_threshold = 1,
          bh_p_value_threshold = 0.05
        )
      })

      output$RESULTS_ORA_TABLE_GO <- DT::renderDT({
        req(
          input$RESULTS_ORA_CATEGORY_CHOICE,
          input$RESULTS_ORA_TYPE_CHOICE
        )
        display_ORA_table(
          ORA_table = .object_@ora_table$GO_TERMS,
          category_choice = "GO_TERMS",
          type_choice = input$RESULTS_ORA_TYPE_CHOICE
        )
      })

      output$RESULTS_ORA_TABLE_KEGG <- DT::renderDT({
        req(
          input$RESULTS_ORA_CATEGORY_CHOICE,
          input$RESULTS_ORA_TYPE_CHOICE
        )
        display_ORA_table(
          ORA_table = .object_@ora_table$KEGG_PWS,
          category_choice = "KEGG_PWS",
          type_choice = input$RESULTS_ORA_TYPE_CHOICE
        )
      })

    })
}

display_parameters_table <- function(
  parameters
) {
  param <- unlist(parameters)
  param <- data.table(
    "Parameter" = names(param),
    "Value" = param
  )
  PARAM_DT <- DT::datatable(
    data = param,
    class = "display compact",
    options =list(
      pageLength = 25,
      dom = '<"top"f>rt<"bottom"lip><"clear">'
    ),
    # caption = tags$caption(
    #   style = paste0(
    #     'caption-side: top; text-align: center; ',
    #     'color:black; font-size:120% ;'
    #   ),
    #   "Caption here"
    # ),
    rownames = rownames,
    extensions = c("Buttons")
  ) #%>% DT::formatStyle(
  #colnames(dt[, -c(9, 10)])[4:8],
  #`text-align` = 'center'
  #)
}

subset_CCI_table <- function(
  CCI_table,
  emitter_choice = NULL,
  receiver_choice = NULL,
  LRI_choice = NULL,
  GENE_choice = NULL,
  GO_choice = NULL,
  KEGG_choice = NULL,
  filter
) {
  dt <- copy(CCI_table)
  if (filter) {
    dt <- dt[
      EMITTER_CELLTYPE %in% emitter_choice &
        RECEIVER_CELLTYPE %in% receiver_choice
    ]
    if (!("All Genes" %in% GENE_choice)) {
      dt <- dt[
        (LIGAND_1 %in% GENE_choice | LIGAND_2 %in% GENE_choice |
           RECEPTOR_1 %in% GENE_choice | RECEPTOR_2 %in% GENE_choice |
           RECEPTOR_3 %in% GENE_choice)
      ]
    }
    if (!("All LRIs" %in% LRI_choice)) {
      dt <- dt[
        LRI %in% LRI_choice
      ]
    }
    if (!("All GO Terms" %in% GO_choice)) {
      GO_ID_choice <- unique(
        .object_@ora_table$GO_TERMS[
          VALUE %in% GO_choice
        ]$VALUE_BIS
      )
      if (.object_@parameters$LRI_species == "human") {
        LRI_GO <- scDiffCom::LRI_human$LRI_curated_GO[
          GO_ID %in% GO_ID_choice
        ]$LRI
      } else if (.object_@parameters$LRI_species == "mouse") {
        LRI_GO <- scDiffCom::LRI_mouse$LRI_curated_GO[
          GO_ID %in% GO_ID_choice
        ]$LRI
      } else if (.object_@parameters$LRI_species == "rat") {
        LRI_GO <- scDiffCom::LRI_rat$LRI_curated_GO[
          GO_ID %in% GO_ID_choice
        ]$LRI
      } else {
        stop(
          "Species",
          .object_@parameters$LRI_species,
          "not supported"
        )
      }
      dt <- dt[
        LRI %in% LRI_GO
      ]
    }
    if (!("All KEGG Pathways" %in% KEGG_choice)) {
      if (.object_@parameters$LRI_species == "human") {
        LRI_KEGG <- scDiffCom::LRI_human$LRI_curated_KEGG[
          KEGG_NAME %in% KEGG_choice
        ]$LRI
      } else if (.object_@parameters$LRI_species == "mouse") {
        LRI_KEGG <- scDiffCom::LRI_mouse$LRI_curated_KEGG[
          KEGG_NAME %in% KEGG_choice
        ]$LRI
      } else if (.object_@parameters$LRI_species == "rat") {
        LRI_KEGG <- scDiffCom::LRI_rat$LRI_curated_KEGG[
          KEGG_NAME %in% KEGG_choice
        ]$LRI
      } else {
        stop(
          "Species",
          .object_@parameters$LRI_species,
          "not supported"
        )
      }
      dt <- dt[
        LRI %in% LRI_KEGG
      ]
    }
  }
  dt
}

plot_volcano_CCI <- function(
  CCI_table,
  params
) {
  dt <- copy(CCI_table)
  dt$Regulation <- factor(
    dt$REGULATION,
    levels = c("UP", "DOWN", "FLAT", "NSC")
  )
  dt$LOG2FC <- dt$LOGFC/log(2)
  dt$MLOG10P <- ifelse(
    dt$BH_P_VALUE_DE == 0,
    -log10(1/params$iterations),
    -log10(dt$BH_P_VALUE_DE)
  )
  vline <- function(x = 0, color = "black") {
    list(
      type = "line",
      y0 = 0,
      y1 = 1,
      yref = "paper",
      x0 = x,
      x1 = x,
      line = list(color = color)
    )
  }
  hline <- function(y = 0, color = "black") {
    list(
      type = "line",
      x0 = 0,
      x1 = 1,
      xref = "paper",
      y0 = y,
      y1 = y,
      line = list(color = color)
    )
  }
  m <- list(
    l = 10,
    r = 10,
    b = 10,
    t = 30,
    pad = 10
  )
  plotly::plot_ly(
    data = dt,
    type = "scatter",
    mode = "markers",
    x = ~LOG2FC,
    y = ~MLOG10P,
    text = ~paste(
      "LRI: ",
      LRI,
      '<br>Emitter:',
      EMITTER_CELLTYPE,
      '<br>Receiver:',
      RECEIVER_CELLTYPE
    ),
    color = ~Regulation,
    colors = stats::setNames(
      c("red", "blue", "green", "gray"),
      c("UP", "DOWN", "FLAT", "NSC")
    )
  ) %>% plotly::layout(
    title = list(
      text = "Interactive Volcano Plot",
      font = list(size = 20),
      xanchor = "left",
      x = 0.0
    ),
    xaxis = list(
      title = list(
        text = "Log2(FC)",
        font = list(size = 18)
      )
    ),
    yaxis = list(
      title = list(
        text = "-Log10(Adj. p-value)",
        font = list(size = 18)
      ),
      range = list(0, -log10(1/params$iterations)*1.2)
    ),
    shapes = list(
      vline(params$threshold_logfc/log(2)),
      vline(-params$threshold_logfc/log(2)),
      hline(-log10(params$threshold_p_value_de))
    ),
    legend = list(
      orientation = "h",
      xanchor = "center",
      x = 0.5,
      y = 1.02
    ),
    margin = m
  )# %>% plotly::toWebGL()
}

plot_scores_CCI <- function(
  CCI_table,
  params
) {
  dt <- copy(CCI_table)
  dt$Regulation <- factor(
    dt$REGULATION,
    levels = c("UP", "DOWN", "FLAT", "NSC")
  )
  cond1_score <- paste0("CCI_SCORE_", params$seurat_condition_id$cond1_name)
  cond2_score <- paste0("CCI_SCORE_", params$seurat_condition_id$cond2_name)
  min_score <-  10^(floor(
    log10(
      min(
        min(dt[get(cond1_score) > 0][[cond1_score]]),
        min(dt[get(cond2_score) > 0][[cond2_score]])
      )
    )
  ))
  dt$x_axis <- ifelse(
    dt[[cond1_score]] == 0,
    log10(min_score),
    log10(dt[[cond1_score]])
  )
  dt$y_axis <- ifelse(
    dt[[cond2_score]] == 0,
    log10(min_score),
    log10(dt[[cond2_score]])
  )
  m <- list(
    l = 10,
    r = 10,
    b = 10,
    t = 30,
    pad = 10
  )
  plotly::plot_ly(
    data = dt,
    type = "scatter",
    mode = "markers",
    x = ~x_axis,
    y = ~y_axis,
    text = ~paste(
      "LRI: ",
      LRI,
      '<br>Emitter:',
      EMITTER_CELLTYPE,
      '<br>Receiver:',
      RECEIVER_CELLTYPE
    ),
    color = ~Regulation,
    colors = stats::setNames(
      c("red", "blue", "green", "gray"),
      c("UP", "DOWN", "FLAT", "NSC")
    )
  ) %>% plotly::layout(
    title = list(
      text = "Interactive Score Plot",
      font = list(size = 20),
      xanchor = "left",
      x = 0.0
    ),
    xaxis = list(
      title = list(
        text = paste0("Log10(", cond1_score, ")"),
        font = list(size = 18)
      )
    ),
    yaxis = list(
      title = list(
        text = paste0("Log10(", cond2_score, ")"),
        font = list(size = 18)
      )
    ),
    legend = list(
      orientation = "h",
      xanchor = "center",
      x = 0.5,
      y = 1.02
    ),
    margin = m
  )# %>% plotly::toWebGL()
}

plot_lrfc_CCI <- function(
  CCI_table,
  params
) {
  dt <- copy(CCI_table)
  dt$Regulation <- factor(
    dt$REGULATION,
    levels = c("UP", "DOWN", "FLAT", "NSC")
  )
  dt[
    ,
    LOG2FC_L := log2(
      do.call(
        pmin,
        c(
          lapply(
            1:params$max_nL,
            function(i) {
              get(
                paste0(
                  "L",
                  i,
                  "_EXPRESSION_",
                  params$seurat_condition_id$cond2_name
                )
              )
            }
          ),
          na.rm = TRUE
        )
      )
      /
        do.call(
          pmin,
          c(
            lapply(
              1:params$max_nL,
              function(i) {
                get(
                  paste0(
                    "L",
                    i,
                    "_EXPRESSION_",
                    params$seurat_condition_id$cond1_name
                  )
                )
              }
            ),
            na.rm = TRUE
          )
        )
    )
  ]
  dt[
    ,
    LOG2FC_R := log2(
      do.call(
        pmin,
        c(
          lapply(
            1:params$max_nR,
            function(i) {
              get(
                paste0(
                  "R",
                  i,
                  "_EXPRESSION_",
                  params$seurat_condition_id$cond2_name
                )
              )
            }
          ),
          na.rm = TRUE
        )
      )
      /
        do.call(
          pmin,
          c(
            lapply(
              1:params$max_nR,
              function(i) {
                get(
                  paste0(
                    "R",
                    i,
                    "_EXPRESSION_",
                    params$seurat_condition_id$cond1_name
                  )
                )
              }
            ),
            na.rm = TRUE
          )
        )
    )
  ]
  dt[
    ,
    LOG2FC_L := {
      max_L <- ceiling(max(.SD[is.finite(LOG2FC_L)][["LOG2FC_L"]]))
      min_L <- floor(min(.SD[is.finite(LOG2FC_L)][["LOG2FC_L"]]))
      max_L <- max(max_L, -min_L)
      min_L <- min(-max_L, min_L)
      ifelse(
        is.infinite(LOG2FC_L) & LOG2FC_L > 0,
        max_L,
        ifelse(
          is.infinite(LOG2FC_L) & LOG2FC_L < 0,
          min_L,
          LOG2FC_L
        )
      )
    }
  ]
  dt[
    ,
    LOG2FC_R := {
      max_R <- ceiling(max(.SD[is.finite(LOG2FC_R)][["LOG2FC_R"]]))
      min_R <- floor(min(.SD[is.finite(LOG2FC_R)][["LOG2FC_R"]]))
      max_R <- max(max_R, -min_R)
      min_R <- min(-max_R, min_R)
      ifelse(
        is.infinite(LOG2FC_R) & LOG2FC_R > 0,
        max_R,
        ifelse(
          is.infinite(LOG2FC_R) & LOG2FC_R < 0,
          min_R,
          LOG2FC_R
        )
      )
    }
  ]
  m <- list(
    l = 10,
    r = 10,
    b = 10,
    t = 30,
    pad = 10
  )
  plotly::plot_ly(
    data = dt,
    type = "scatter",
    mode = "markers",
    x = ~LOG2FC_L,
    y = ~LOG2FC_R,
    text = ~paste(
      "LRI: ",
      LRI,
      '<br>Emitter:',
      EMITTER_CELLTYPE,
      '<br>Receiver:',
      RECEIVER_CELLTYPE
    ),
    color = ~Regulation,
    colors = stats::setNames(
      c("red", "blue", "green", "gray"),
      c("UP", "DOWN", "FLAT", "NSC")
    )
  )  %>% plotly::layout(
    title = list(
      text = "Interactive 'Ligand-FC vs Receptor-FC' Plot",
      font = list(size = 20),
      xanchor = "left",
      x = 0.0
    ),
    xaxis = list(
      title = list(
        text = "Ligand Log2(FC)",
        font = list(size = 18)
      )
    ),
    yaxis = list(
      title = list(
        text = "Receptor Log2(FC)",
        font = list(size = 18)
      )
    ),
    legend = list(
      orientation = "h",
      xanchor = "center",
      x = 0.5,
      y = 1.02
    ),
    margin = m
  )# %>% plotly::toWebGL()
}

display_CCI_table <- function(
  CCI_table,
  params
) {
  dt <- copy(CCI_table)
  dt[, LOG2FC := signif(dt$LOGFC/log(2), 3)]
  dt[, BH_P_VALUE_DE := signif(BH_P_VALUE_DE, 3)]
  dt[
    ,
    paste0(
      "CCI_SCORE_",
      params$seurat_condition_id$cond1_name
    ) :=
      signif(
        get(
          paste0(
            "CCI_SCORE_",
            params$seurat_condition_id$cond1_name
          )
        ),
        4
      )
  ]
  dt[
    ,
    paste0(
      "CCI_SCORE_",
      params$seurat_condition_id$cond2_name
    ) :=
      signif(
        get(
          paste0(
            "CCI_SCORE_",
            params$seurat_condition_id$cond2_name
          )
        ),
        4
      )
  ]
  cols <- c(
    "EMITTER_CELLTYPE",
    "RECEIVER_CELLTYPE",
    "LRI",
    paste0("CCI_SCORE_", params$seurat_condition_id$cond1_name),
    paste0("CCI_SCORE_", params$seurat_condition_id$cond2_name),
    "LOG2FC",
    "BH_P_VALUE_DE",
    "REGULATION"
  )
  dt <- dt[
    ,
    ..cols
  ]
  CCI_DT <- DT::datatable(
    data = dt,
    class = "display compact",
    options =list(
      pageLength = 10,
      dom = '<"top"f>rt<"bottom"lip><"clear">'
    ),
    caption = tags$caption(
      style = paste0(
        'caption-side: top; text-align: center; ',
        'color:black; font-size:120% ;'
      ),
      "Table of Cell-Cell Interactions"
    ),
    rownames = rownames,
    extensions = c("Buttons")
  ) %>% DT::formatStyle(
    cols,
    `text-align` = 'center'
  )
}

display_ORA_table <- function(
  ORA_table,
  category_choice,
  type_choice
) {
  dt <- copy(ORA_table)
  if (category_choice == "GO_TERMS") {
    level_str <- c("LEVEL", "ASPECT")
    filter <- "top"
  } else {
    level_str <- NULL
    filter <- "none"
  }
  if(type_choice == "UP") {
    cols_to_keep <- c(
      "VALUE",
      "ORA_SCORE_UP",
      "OR_UP",
      "BH_P_VALUE_UP",
      level_str
    )
    dt <- dt[
      OR_UP >= 1 & BH_P_VALUE_UP <= 0.05,
      cols_to_keep,
      with = FALSE
    ]
  } else if(type_choice == "DOWN") {
    cols_to_keep <- c(
      "VALUE",
      "ORA_SCORE_DOWN",
      "OR_DOWN",
      "BH_P_VALUE_DOWN",
      level_str
    )
    dt <- dt[
      `OR_DOWN` >= 1 & BH_P_VALUE_DOWN <= 0.05,
      cols_to_keep,
      with = FALSE
    ]
  } else if(type_choice == "FLAT") {
    cols_to_keep <- c(
      "VALUE",
      "ORA_SCORE_FLAT",
      "OR_FLAT",
      "BH_P_VALUE_FLAT",
      level_str
    )
    dt <- dt[
      `OR_FLAT` >= 1 & BH_P_VALUE_FLAT <= 0.05,
      cols_to_keep,
      with = FALSE]
  }
  if (category_choice == "GO_TERMS") {
    dt[, LEVEL := as.factor(LEVEL)]
    dt[, ASPECT := as.factor(ASPECT)]
  }
  data.table::setnames(
    dt,
    old = colnames(dt),
    new = c(
      category_choice,
      "ORA Score",
      "Odds Ratio",
      "Adj. p-value",
      level_str
    )
  )
  data.table::setorder(dt, -`ORA Score`)
  dt[, `ORA Score` := signif(`ORA Score`, 4)]
  dt[, `Odds Ratio` := signif(`Odds Ratio`, 4)]
  dt[, `Adj. p-value` := signif(`Adj. p-value`, 3)]
  DT::datatable(
    data = dt,
    options = list(
      columnDefs = list(
        list(width = '300px', targets = c(0))
      ),
      dom = '<"top"f>rt<"bottom"lip><"clear">'
    ),
    caption = tags$caption(
      style = paste0(
        "caption-side: top; ",
        "text-align: center; ",
        "color: black; ",
        "font-size: 120%;"
      ),
      paste0(
        category_choice,
        " over-represented among ",
        type_choice,
        "-regulated cell-cell interactions"
      )
    ),
    rownames = FALSE,
    filter = filter,
    class = "display compact"
  ) %>% DT::formatStyle(
    colnames(dt)[-1],
    `text-align` = 'center'
  )
}

plot_ORA_GO_treemap <- function(
  GO_REDUCED_table,
  type_choice,
  go_aspect_choice,
  title_text,
  domain = NULL
) {
  ex_data <- GO_REDUCED_table[
    ASPECT == go_aspect_choice &
      REGULATION == type_choice
  ][, c("score", "term", "parentTerm")]
  if (nrow(ex_data) == 0) return(NULL)
  ex_data[, new_parent := ifelse(
    term %in% parentTerm,
    "",
    parentTerm
  )]
  new_data <- data.table(
    labels = c(ex_data$term, ex_data[new_parent == ""]$term),
    parents = c(
      ex_data$parentTerm,
      rep("", length(ex_data[new_parent == ""]$term))
    )
  )
  new_data[
    ,
    ids := sapply(
      1:nrow(.SD),
      function(i) {
        if (labels[[i]] == parents[[i]]) {
          res <- paste(labels[[i]], parents[[i]], sep = " - ")
        } else {
          res <- labels[[i]]
        }
        res
      }
    )
  ]
  new_data[
    ,
    score := sapply(
      1:nrow(.SD),
      function(i) {
        if (parents[[i]] == "") {
          res <- sum(ex_data[parentTerm == labels[[i]]]$score)
        } else {
          res <- ex_data[term == labels[[i]]]$score
        }
        res
      }
    )
  ]
  new_data[
    ,
    text := gsub(" ", "\n", labels)
  ]
  m <- list(
    l = 5,
    r = 5,
    b = 5,
    t = 30,
    pad = 0
  )
  plotly::plot_ly(
    new_data,
    type = "treemap",
    opacity = 1,
    ids = ~ids,
    parents = ~parents,
    values = ~score,
    labels = ~labels,
    text = ~text,
    textposition = "middle center",
    branchvalues = "total",
    hoverinfo = "label+value",
    marker = list(
      line = list(color = "black")
    ),
    textinfo = "text",
    domain = domain
  ) %>% plotly::layout(
    title = list(
      text = title_text,
      font = list(size = 16),
      xanchor = "left",
      x = 0.0
    ),
    #uniformtext = list(
    #  minsize = 14,
    #  mode = "hide"
    #),
    margin = m
  )
}
