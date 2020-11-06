#
# build_celltype_bipartite_graph <- function(
#   object
# ) {
#   ORA_tables <- get_ora_tables(object)
#   if("LR_CELLTYPE" %in% names(ORA_tables)) {
#     ORA_ct <- ORA_tables[["LR_CELLTYPE"]]
#   } else {
#     object <- run_ORA(
#       object = object,
#       categories = "LR_CELLTYPE"
#     )
#     ORA_ct <- get_ora_tables(object)[["LR_CELLTYPE"]]
#   }
#   graph <- construct_graph(
#     ORA_ct = ORA_ct
#     )
#   config <- setup_graph_config()
#   graph <- setup_graph(
#     graph = graph,
#     config = config,
#     cci_table_filtered = get_cci_table_filtered(object)
#   )
#
#
#   return(graph)
#
#
#
#
#
#   #if( is.null(dt_filtered) ) {stop("analyze_Graph: dt_filtered is NULL.")}
#   # Process ORA results and construct graph
#   if ( !is.null(dir) ) {
#     if ( is.null(analysis_name) ) {
#       stop("analyze_Graph: Analysis name not provided.")
#     }
#     # save_run_metadata()
#     subdirs = c("edge_tables", "plots")
#     create_analysis_dirs(dir, analysis_name, subdirs)
#     write_as_edge_table(
#       graph,
#       path = file.path(dir, analysis_name, "edge_tables", paste0(tissue, ".tsv"))
#     )
#     plot_graph(
#       graph,
#       config,
#       path = file.path(dir, analysis_name, "plots", paste0(tissue, ".pdf"))
#     )
#   } else {
#     if ( !is.null(analysis_name) ) {
#       stop("analyze_Graph: Analysis name not null.")
#     }
#     plot_graph(
#       graph,
#       config,
#       path=NULL)
#   }
# }
#
# construct_graph <- function(
#   ORA_ct
# ) {
#   ORA_types <- c("UP", "DOWN", "DIFF")
#   dt_edge <- get_celltypes_enrichment(
#     ORA_ct = ORA_ct,
#     ORA_types = ORA_types
#   )
#   dt_edge <- dt_edge[, list(
#     "L_CELLTYPE" = paste0(L_CELLTYPE, " (L)"),
#     "R_CELLTYPE" = paste0(R_CELLTYPE, " (R)")
#   )]
#   graph <- igraph::graph_from_data_frame(
#     d = dt_edge,
#     directed = TRUE,
#     vertices = NULL
#   )
#   return(graph)
#   #graph$interacts_counts <- cci_table_filtered[
#   #  list(num_interacts = .N),
#   #  by = .(L_CELLTYPE, R_CELLTYPE)
#   #  ]
#   #graph$ora <- ORA_ct
#   #graph$dt <- cci_table_filtered
#   #
#   # dt_edge <- dt_edge[, list(
#   #   "L_CELLTYPE" = paste0(L_CELLTYPE, " (L)"),
#   #   "R_CELLTYPE" = paste0(R_CELLTYPE, " (R)"),
#   #   "OR_DIFF" = OR_DIFF,
#   #   "OR_UP" = OR_UP,
#   #   "OR_DOWN" = OR_DOWN,
#   #   "pval_DIFF" = pval_DIFF,
#   #   "pval_UP" = pval_UP,
#   #   "pval_DOWN" = pval_DOWN,
#   #   "pval_adj" = pval_readj_on_celltypes,
#   #   "pval_adj_UP" =  pval_readj_on_celltypes_UP,
#   #   "pval_adj_DOWN" = pval_readj_on_celltypes_DOWN
#   # )]
#
#   # Handle NAs
#   #message(paste0("process_edge_dt: #NAs in ORA:celltypes for ",
#   #               tissue, ":", sum(is.na(dt_edge)),"; NA->1"))
#   #dt_edge[is.na(dt_edge)] = 1
#   # Handle Inf
#   # Shouldn"t have Inf with pseudocounts
#   # CLIP_VAL = 50
#   #if (sum(dt_edge==Inf) > 0) {
#   #  stop("construct_graph: Inf values in dt_edge")
#     # message(paste0("process_edge_dt: Inf value present in edge data.table for",
#     #                tissue, "; Clipping to ", CLIP_VAL))
#     # dt_edge[dt_edge == Inf] = min(dt_edge[, list(OR, OR_UP, OR_DOWN)], CLIP_VAL)
#   #}
# }
#
# setup_graph_config <- function(
# ) {
#   GRAPH_CONFIG = list(
#     EDGE_COLORING = list(
#       COLOR_UP = "#4285F4",
#       CUTOFF_UP = 1.1,
#       COLOR_DOWN = "#DB4437",
#       CUTOFF_DOWN = 1.1,
#       COLOR_BOTH = "#F4B400",
#       COLOR_ROBUST = "#0F9D58",
#       # COLOR_NONE = "#dde2e4"
#       COLOR_NONE = rgb(0.2, 0.2, 0.2, alpha=0.2)
#     ),
#     EDGE_STYLE = list(
#       ARROW_SIZE = 0.5,
#       WIDTH = 2.5  # will be adjusted by |OR|
#     ),
#     LAYOUT = list(
#       HGAP = 10,
#       VGAP = 10
#     ),
#     VERTEX_STYLE = list(
#       SIZE = 10,
#       LABEL_DIST = 1.5,
#       LABEL_CEX = 1.2,
#       COLOR = '#33FF66'
#     ),
#     LEGEND = list(
#       LEGEND_LABELS = c(
#         'Significant, but small effect',
#         'Upregulated',
#         'Downregulated',
#         'Altered'),
#       PCH = c(15),
#       CEX = 0.7,
#       PT.CEX = 1,
#       BG = '#CCCCCC',
#       NCOL = 2
#     )
#   )
# }
#
# setup_graph <- function(
#   graph,
#   config,
#   cci_table_filtered
# ) {
#   igraph::V(graph)$vertex_types <- grepl("(L)", igraph::V(graph)$name, fixed = TRUE)
#   graph <- add_vertex_size(
#     graph,
#     cci_table_filtered
#     )
#   igraph::E(graph)$arrow.size <- config$EDGE_STYLE$ARROW_SIZE
#   igraph::E(graph)$edge.color <- purrr::map_chr(igraph::E(graph), ~ color_edge(.x, config))
#   return(graph)
#
#
#   graph = add_edge_width(graph)
#   layout = layout_as_bipartite(
#     graph,
#     types = V(graph)$vertex_types,
#     hgap = config$LAYOUT$HGAP,
#     vgap = config$LAYOUT$VGAP
#   )
#   layout = layout[, 2:1]  # horizontal to vertical
#   graph$layout = layout
#   return(graph)
# }
#
# add_vertex_size <- function(
#   graph,
#   cci_table_filtered
# ) {
#   igraph::V(graph)$vertex.size <- 30
#   MAXSIZE <- 20
#   MINSIZE <- 5
#   vertex_names <- igraph::V(graph)$name
#   celltypes_to_vertexnames <- function(counts_dt) {
#     vertex_dt = funion(
#       counts_dt[, .(
#         name = paste0(name, ' (L)'),
#         counts = counts
#       )],
#       counts_dt[, .(
#         name = paste0(name, ' (R)'),
#         counts = counts
#       )]
#     )
#     return(vertex_dt)
#   }
#   graph$celltype_counts <- count_celltypes(cci_table_filtered)
#   vertex_counts_dt <- celltypes_to_vertexnames(graph$celltype_counts)
#   vertex_counts_dt <-  vertex_counts_dt[vertex_counts_dt$name %in% vertex_names]
#   ord <- match(vertex_names, vertex_counts_dt$name)
#   vertex_counts <- vertex_counts_dt[ord, counts]
#   igraph::V(graph)$vertex.counts = vertex_counts
#   vertex_sizes = scales::rescale(
#     igraph::V(graph)$vertex.counts,
#     to = c(MINSIZE, MAXSIZE)
#   )
#   igraph::V(graph)$vertex.size = vertex_sizes
#   return(graph)
# }
#
# count_celltypes <- function(
#   cci_table_filtered
# ) {
#   dt_t <- data.table::copy(cci_table_filtered)
#   counts_dt <- data.table::funion(
#     dt_t[, .(
#       name = L_CELLTYPE,
#       counts = (L_NCELLS_OLD + L_NCELLS_YOUNG)/2  # average counts in young-old
#     )],
#     dt_t[, .(
#       name = R_CELLTYPE,
#       counts = (R_NCELLS_OLD + R_NCELLS_YOUNG)/2  # average counts in young-old
#     )]
#   )
#   counts_dt <- counts_dt[order(counts_dt$name)]
#   return(counts_dt)
# }
#
# color_edge <- function(
#   edge,
#   config
# ) {
#   COLOR_UP_AND_DOWN = config$EDGE_COLORING$COLOR_BOTH
#   COLOR_UP = config$EDGE_COLORING$COLOR_UP
#   COLOR_DOWN = config$EDGE_COLORING$COLOR_DOWN
#   COLOR_ROBUST = config$EDGE_COLORING$COLOR_ROBUST
#   COLOR_NONE = config$EDGE_COLORING$COLOR_NONE
#
#   CUTOFF_CHANGE_UP = config$EDGE_COLORING$CUTOFF_UP
#   CUTOFF_CHANGE_DOWN = config$EDGE_COLORING$CUTOFF_DOWN
#   CUTOFF_ROBUST_UP = 0.99
#   CUTOFF_ROBUST_DOWN = 0.99
#   #message('color_edge: some hard-coded cutoff.')
#   #message(paste0(
#   #  'color_edge: refactor to set edge state in the ',
#   #  'edge data table to be reused within further analysis.')
#   #)
#   #if (use_adjpval) {
#   is_pval_UP = edge$pval_adj_UP < 0.05
#   is_pval_DOWN = edge$pval_adj_DOWN < 0.05
#   #} else {
#   #  is_pval_UP = edge$pval_UP < 0.05
#   #  is_pval_DOWN = edge$pval_DOWN < 0.05
#   #}
#   is_OR_UP_CHANGE = edge$OR_UP >= CUTOFF_CHANGE_UP
#   is_OR_DOWN_CHANGE = edge$OR_DOWN >= CUTOFF_CHANGE_DOWN
#   is_OR_UP_ROBUST = edge$OR_UP <= CUTOFF_ROBUST_UP
#   is_OR_DOWN_ROBUST = edge$OR_DOWN <= CUTOFF_ROBUST_DOWN
#     if (is_OR_UP_CHANGE & is_OR_DOWN_CHANGE & is_pval_UP & is_pval_DOWN) {
#     return(COLOR_UP_AND_DOWN)
#   } else if (is_OR_UP_CHANGE & is_pval_UP) {
#     return(COLOR_UP)
#   } else if (is_OR_DOWN_CHANGE & is_pval_DOWN) {
#     return(COLOR_DOWN)
#   } else if (is_OR_UP_ROBUST & is_OR_DOWN_ROBUST & is_pval_UP & is_pval_DOWN) {
#     return(COLOR_ROBUST)
#   } else {
#     return(COLOR_NONE)
#   }
# }
#
# ############
#
#
#
# add_edge_width <- function(graph) {
#
#   WIDTH_MIN = 0
#   WIDTH_MAX = 10
#
#   # or_max = mapply(max, E(graph)$OR_UP, E(graph)$OR_DOWN)
#
#   edge_dt = data.table(
#     as_data_frame(graph)[, c('from', 'to')]
#   )
#   names(edge_dt) = c('Ligand_celltype', 'Receptor_celltype')
#
#   edge_dt[, Ligand_celltype :=
#             sub('(*) [(]L[)]', '', Ligand_celltype)
#           ][,
#             Receptor_celltype :=
#               sub('(*) [(]R[)]', '', Receptor_celltype)
#             ]
#
#   edge_dt = merge(
#     edge_dt, graph$interacts_counts,
#     by = c('Ligand_celltype', 'Receptor_celltype'),
#     all.x = TRUE,
#     sort = FALSE)
#
#   edge_dt[is.na(edge_dt)] = 0
#   message('add_edge_width: Some NAs occur, unclear reason. For now NA->0.')
#
#   # Correct order is accomplished by the merge.
#   num_interacts = edge_dt[, num_interacts]
#
#   E(graph)$width = rescale(
#     sqrt(num_interacts),
#     to = c(WIDTH_MIN, WIDTH_MAX)
#   )
#   message('add_edge_width: widths are scaled non-linearly.')
#
#   return(graph)
# }
#
#
# get_celltypes_enrichment <- function(
#   ORA_ct,
#   ORA_types
# ) {
#   cols_pval <- paste0("pval_", ORA_types)
#   cols_pval_readj <- paste0("pval_readj_on_celltypes_", ORA_types)
#   dt_ctypes <- ORA_ct[
#     Category == "LR_CELLTYPE"
#     ][,
#       (cols_pval_readj) := lapply(
#         .SD,
#         stats::p.adjust,
#         method = "BH"
#       ),
#       .SDcols = cols_pval
#       ]
#   dt_ctypes <- dt_ctypes[, c("L_CELLTYPE", "R_CELLTYPE") := list(
#     sub("_.*", "", Value),  # substitute w/ regex
#     sub(".*_", "", Value)
#   )
#   ]
#   cols_to_select = c(
#     "L_CELLTYPE", "R_CELLTYPE", paste0("OR_", ORA_types),
#     cols_pval, cols_pval_readj
#   )
#   dt_ctypes <- dt_ctypes[, cols_to_select, with = FALSE]
#   return(dt_ctypes)
# }
#
# plot_graph <- function(G,
#                        config=GRAPH_CONFIG,
#                        path=NULL,
#                        show_legend=FALSE) {
#
#   # Plot
#   MAIN = paste0('Communication profile with age: ', G$tissue)
#   MARGIN = 0.5
#   SUBTITLE = ''
#   # SUBTITLE = paste0('Difference graph of overrepresented celltype ',
#   #                   'communication that get\n altered with ageing')
#
#   if( !is.null(path) ) {pdf(path)}
#
#   plot.igraph(
#     G,
#     layout = G$layout,
#     vertex.color = config$VERTEX_STYLE$COLOR,
#     vertex.size = V(G)$vertex.size,
#     vertex.label.dist = config$VERTEX_STYLE$LABEL_DIST,
#     vertex.label.cex = config$VERTEX_STYLE$LABEL_CEX,
#     edge.color = E(G)$edge.color,
#     edge.width = E(G)$edge.width,
#     main = MAIN,
#     margin = MARGIN,
#     sub = SUBTITLE
#   )
#
#   if (show_legend) {
#     LEGEND_TITLE = 'Edge color legend'
#     LEGEND_COLORS = c(
#       config$EDGE_COLORING$COLOR_NONE,
#       config$EDGE_COLORING$COLOR_UP,
#       config$EDGE_COLORING$COLOR_DOWN,
#       config$EDGE_COLORING$COLOR_BOTH
#     )
#
#     legend(x=-1.5, y=-1.1,
#
#            legend = config$LEGEND$LEGEND_LABELS,
#            col = LEGEND_COLORS,
#            title = LEGEND_TITLE,
#
#            bty = 'o',
#            pch = config$LEGEND$PCH,
#            cex = config$LEGEND$CEX,
#            pt.cex = config$LEGEND$PT.CEX,
#            bg = config$LEGEND$BG,
#            ncol = config$LEGEND$NCOL
#     )
#   }
#   if( !is.null(path) ) {dev.off()}
# }
#
# #####################
#
# build_celltype_heatmaps <- function(
#   ORA_dt,
#   ORA_types
#   ) {
#   OR_max <- 10
#   OR_min <- 0.1
#   #eps <- 1.E-5
#   dt_edge <- get_celltypes_enrichment(
#     ORA_dt = ORA_dt,
#     ORA_types = ORA_types
#   )
#   graph <-  igraph::graph_from_data_frame(
#     d = dt_edge,
#     directed = TRUE,
#     vertices = NULL
#   )
#   adja_matrix <- lapply(
#     ORA_types,
#     function(type) {
#       temp <- igraph::as_adjacency_matrix(
#         graph = graph,
#         attr = paste0("OR_", type),
#         sparse = FALSE
#       )
#       temp <- clip_matrix(
#         matrix = temp,
#         zero_val = 1,
#         min_val = OR_min,
#         max_val = OR_max
#       )
#     }
#   )
#   plots_hm <- lapply(
#     adja_matrix,
#     function(mat) {
#       generate_heatmap(
#         matrix = log(mat),
#         #matrix = log(mat + eps),
#         title = "test"
#       )
#
#     }
#   )
#
#
#
#   #nrows = dim(dt_edge)[1]
#   #if (nrows == 1) {message(paste0("One edge for: ", tissue))}
#   #dir.create(paste0(dir, tissue, "/"), showWarnings = FALSE)
#
#   #explanatory_string = paste0("Row-transmitter, col-receiver; ",
#   #                            "Values are clipped log(OR)")
#
#   #generate_heatmap(log(matrix_all + eps),
#   #                 #matrix_all,
#   #                 title = paste0(tissue, ": ", "ALL ", explanatory_string),
#   #                 filename = paste0(dir, tissue, "/", "ALL.png"))
#
# }
#
# generate_heatmap <- function(
#   matrix,
#   title#,
#   #filename
# ) {
#   nrows <- dim(matrix)[1]
#   ncols <- dim(matrix)[2]
#   if (nrows < 2 | ncols < 2) {
#     message(paste0("Matrix size is < 2x2: ", title))
#     return()
#   }
#   breaks <- c(seq(-3, -0.5, 0.5), 0, seq(0.5, 3, 0.5))
#   temp_color <- grDevices::colorRampPalette(
#     RColorBrewer::brewer.pal(8, "RdYlBu"))(length(breaks)-1)
#   hm <- pheatmap::pheatmap(
#     mat = matrix,
#     #cellwidth=10,
#     #cellheight=10,
#     width = 11,
#     height = 11,
#     color = temp_color,
#     breaks = breaks,
#     na_col = "green",
#     scale = "none",
#     cluster_rows = nrows >= 2,
#     cluster_cols = ncols >= 2,
#     cutree_rows = ifelse(nrows >= 5, min(nrows, 3), 1),
#     cutree_cols = ifelse(ncols >= 5, min(ncols, 3), 1),
#     display_numbers = FALSE,
#     fontsize = 10,
#     main = title,
#     ylab = "Emitter Cell Type",
#     xlab = "Receiver Cell Type",
#     legend = TRUE#,
#     #filename = filename
#   )
#   return(hm)
# }
#
#
#
# clip_matrix <- function(
#   matrix,
#   zero_val,
#   min_val,
#   max_val
# ) {
#   matrix[matrix == 0] <- zero_val
#   matrix[matrix < min_val] <- min_val
#   matrix[matrix > max_val] <- max_val
#   return(matrix)
# }
#
