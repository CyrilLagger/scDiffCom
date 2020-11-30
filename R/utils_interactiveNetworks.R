# get_celltypes_enrichment <- function(dt) {
#
#     COL_LIGAND_RECEPTOR_CELLTYPES = "ER_CELLTYPES"
#     COL_pval = 'P_VALUE_DIFF'
#     COL_P_VALUE_UP = 'P_VALUE_UP'
#     COL_P_VALUE_DOWN = 'P_VALUE_DOWN'
#
#     # These indicate p-values that have been readjusted
#     #  at celltype level for keeping the signal.
#     COL_pval_readj = 'pval_readj_on_celltypes'
#     COL_pval_readj_UP = 'pval_readj_on_celltypes_UP'
#     COL_pval_readj_DOWN = 'pval_readj_on_celltypes_DOWN'
#
#     dt_ctypes = dt[
#         CATEGORY == COL_LIGAND_RECEPTOR_CELLTYPES
#         ][,
#           c(COL_pval_readj, COL_pval_readj_UP, COL_pval_readj_DOWN) := .(
#               stats::p.adjust(get(COL_pval), method='BH'),
#               stats::p.adjust(get(COL_P_VALUE_UP), method='BH'),
#               stats::p.adjust(get(COL_P_VALUE_DOWN), method='BH')
#           )
#           ]
#
#     dt_ctypes = dt_ctypes[, c("Ligand_cell", "Receptor_cell") := list(
#         sub("_.*", "", VALUE),  # substitute w/ regex
#         sub(".*_", "", VALUE)
#     )
#     ]
#
#     cols_to_select = c(
#         "Tissue", "Ligand_cell", "Receptor_cell", "OR_DIFF", "OR_UP", "OR_DOWN",
#         COL_pval, COL_P_VALUE_UP, COL_P_VALUE_DOWN,
#         COL_pval_readj, COL_pval_readj_UP, COL_pval_readj_DOWN
#     )
#
#     dt_ctypes = dt_ctypes[, ..cols_to_select]
#
#     return(dt_ctypes)
# }
# classify_edges <- function(dt_edge,
#                            config=setup_graph_config(),
#                            use_adjpval=TRUE) {
#
#     LABEL_DIFF = 'diff'  # label UP or DOWN
#     LABEL_UP = 'up'
#     LABEL_DOWN = 'down'
#     LABEL_ROBUST = 'robust'
#     LABEL_NONE = 'none'
#
#     CUTOFF_CHANGE_UP = config$EDGE_COLORING$CUTOFF_UP
#     CUTOFF_CHANGE_DOWN = config$EDGE_COLORING$CUTOFF_DOWN
#     CUTOFF_ROBUST_UP = 0.99
#     CUTOFF_ROBUST_DOWN = 0.99
#
#     edges_classification = purrr::pmap_chr(
#         dt_edge,
#         function(OR_DIFF, OR_UP, OR_DOWN,
#                  P_VALUE_DIFF, P_VALUE_UP, P_VALUE_DOWN,
#                  pval_adj, pval_adj_UP, pval_adj_DOWN,
#                  ...) {
#
#             if (use_adjpval) {
#                 is_P_VALUE_UP = pval_adj_UP < 0.05
#                 is_P_VALUE_DOWN = pval_adj_DOWN < 0.05
#             } else {
#                 is_P_VALUE_UP = P_VALUE_UP < 0.05
#                 is_P_VALUE_DOWN = P_VALUE_DOWN < 0.05
#             }
#
#             is_OR_UP_CHANGE = OR_UP >= CUTOFF_CHANGE_UP
#             is_OR_DOWN_CHANGE = OR_DOWN >= CUTOFF_CHANGE_DOWN
#             is_OR_UP_ROBUST = OR_UP <= CUTOFF_ROBUST_UP
#             is_OR_DOWN_ROBUST = OR_DOWN <= CUTOFF_ROBUST_DOWN
#
#             if (is_OR_UP_CHANGE & is_OR_DOWN_CHANGE & is_P_VALUE_UP & is_P_VALUE_DOWN) {
#                 return(LABEL_DIFF)
#             } else if (is_OR_UP_CHANGE & is_P_VALUE_UP) {
#                 return(LABEL_UP)
#             } else if (is_OR_DOWN_CHANGE & is_P_VALUE_DOWN) {
#                 return(LABEL_DOWN)
#             } else if (is_OR_UP_ROBUST & is_OR_DOWN_ROBUST & is_P_VALUE_UP & is_P_VALUE_DOWN) {
#                 return(LABEL_ROBUST)
#             } else {
#                 return(LABEL_NONE)
#             }
#         }
#     )
# }
# map_edgetype_to_color = function(edge_label, config=setup_graph_config()) {
#     COLOR_DIFF = config$EDGE_COLORING$COLOR_BOTH
#     COLOR_UP = config$EDGE_COLORING$COLOR_UP
#     COLOR_DOWN = config$EDGE_COLORING$COLOR_DOWN
#     COLOR_ROBUST = config$EDGE_COLORING$COLOR_ROBUST
#     COLOR_NONE = config$EDGE_COLORING$COLOR_NONE
#
#     label_to_color_map = list(
#         'diff'=COLOR_DIFF,
#         'up'=COLOR_UP,
#         'down'=COLOR_DOWN,
#         'robust'=COLOR_ROBUST,
#         'none'=COLOR_NONE
#     )
#     return(label_to_color_map[[edge_label]])
# }
# construct_graph <- function(dt_ora, tissue, dt_filtered=NULL) {
#
#     dt_ctypes = get_celltypes_enrichment(dt_ora)
#     dt_edge = dt_ctypes[Tissue == tissue, .SD, .SDcols = !c('Tissue')]
#
#     # Label the transmitter (L) and receiver (R) cells
#     dt_edge = dt_edge[, .(
#         'Ligand_cell' = paste0(Ligand_cell, ' (L)'),
#         'Receptor_cell' = paste0(Receptor_cell, ' (R)'),
#         # 'OR' = OR,
#         'OR_DIFF' = OR_DIFF,
#         'OR_UP' = OR_UP,
#         'OR_DOWN' = OR_DOWN,
#         # 'pval' = pval,
#         'P_VALUE_DIFF' = P_VALUE_DIFF,
#         'P_VALUE_UP' = P_VALUE_UP,
#         'P_VALUE_DOWN' = P_VALUE_DOWN,
#         'pval_adj' = pval_readj_on_celltypes,
#         'pval_adj_UP' =  pval_readj_on_celltypes_UP,
#         'pval_adj_DOWN' = pval_readj_on_celltypes_DOWN
#     )]
#
#     # Handle NAs
#     message(paste0('process_edge_dt: #NAs in ORA:celltypes for ',
#                    tissue, ':', sum(is.na(dt_edge)),'; NA->1'))
#     dt_edge[is.na(dt_edge)] = 1
#
#     # Shouldn't have Inf with pseudocounts
#     if (sum(dt_edge==Inf) > 0) {
#         stop('construct_graph: Inf values in dt_edge')
#     }
#
#     # Annotate edges
#     dt_edge[, edge_type := classify_edges(dt_edge)]
#
#     G = igraph::graph_from_data_frame(d = dt_edge,
#                                       directed = TRUE,
#                                       vertices = NULL)
#
#     # Set graph tissue attribute
#     G$tissue = tissue
#     G$ora = dt_ora
#     G$dt = dt_filtered
#     G$interacts_counts = count_interactions_cellpairs_tissue(dt_filtered, tissue)
#     G$celltype_counts = count_celltypes_tissue(dt_filtered, tissue)
#
#     message('construct_graph: igraph G holds as graph attributes dt_ora and dt_filtered.')
#     return(G)
# }
# infer_vertex_types <- function(vertex_names) {
#     vertex_types = purrr::map_lgl(
#         strsplit(vertex_names, '[()]'),
#         ~ .x[2] == 'L'
#     )
#     return(vertex_types)
# }
# setup_graph <- function(G, config=setup_graph_config(), use_adjpval, disperse) {
#     G = setup_vertices(G, config)
#     G = setup_edges(G, config, use_adjpval)
#     G = setup_layout(G, config, disperse)
#     return(G)
# }
# count_celltypes_tissue <- function(dt, tissue) {
#     # @dt - dt_filtered
#
#     # dt_t = dt[TISSUE == tissue]
#     dt_t = data.table::copy(dt)
#
#     counts_dt = data.table::funion(
#         dt_t[, .(
#             name = EMITTER_CELLTYPE,
#             counts = (EMITTER_NCELLS_OLD + EMITTER_NCELLS_YOUNG)/2  # average counts in young-old
#         )],
#         dt_t[, .(
#             name = RECEIVER_CELLTYPE,
#             counts = (RECEIVER_NCELLS_OLD + RECEIVER_NCELLS_YOUNG)/2  # average counts in young-old
#         )],
#     )
#     counts_dt = counts_dt[order(counts_dt$name)]
#     return(counts_dt)
# }
# add_vertex_size <- function(G) {
#
#     igraph::V(G)$vertex.size = 30
#     MAXSIZE = 20
#     MINSIZE = 5
#
#     vertex_names = igraph::V(G)$name
#
#     celltypes_to_vertexnames <- function(counts_dt) {
#
#         vertex_dt = data.table::funion(
#             counts_dt[, .(
#                 name = paste0(name, ' (L)'),
#                 counts = counts
#             )],
#             counts_dt[, .(
#                 name = paste0(name, ' (R)'),
#                 counts = counts
#             )]
#         )
#
#         return(vertex_dt)
#     }
#
#     vertex_counts_dt = celltypes_to_vertexnames(G$celltype_counts)
#     vertex_counts_dt = vertex_counts_dt[vertex_counts_dt$name %in% vertex_names]
#
#     # Order data.table by vertex_names order and extract counts in correct order.
#     ord = match(vertex_names, vertex_counts_dt$name)
#     vertex_counts = vertex_counts_dt[ord, counts]
#
#     igraph::V(G)$vertex.counts = vertex_counts
#     vertex_sizes = rescale(vertex_counts, MINSIZE, MAXSIZE)
#
#     igraph::V(G)$vertex.size = vertex_sizes
#
#     return(G)
# }
# count_interactions_cellpairs_tissue <- function(dt_filtered, tissue) {
#     return(
#         dt_filtered[,
#                     .(count = .N),
#                     by = .(EMITTER_CELLTYPE, RECEIVER_CELLTYPE)
#                     ][, .(
#                         'Ligand_celltype' = EMITTER_CELLTYPE,
#                         'Receptor_celltype' = RECEIVER_CELLTYPE,
#                         'num_interacts' = count
#                     )]
#     )
# }
# add_edge_width <- function(G) {
#
#     WIDTH_MIN = 0
#     WIDTH_MAX = 10
#
#     edge_dt = data.table::data.table(
#         igraph::as_data_frame(G)[, c('from', 'to')]
#     )
#     names(edge_dt) = c('Ligand_celltype', 'Receptor_celltype')
#
#     edge_dt[, Ligand_celltype :=
#                 sub('(*) [(]L[)]', '', Ligand_celltype)
#             ][,
#               Receptor_celltype :=
#                   sub('(*) [(]R[)]', '', Receptor_celltype)
#               ]
#
#     edge_dt = merge(
#         edge_dt, G$interacts_counts,
#         by = c('Ligand_celltype', 'Receptor_celltype'),
#         all.x = TRUE,
#         sort = FALSE)
#
#     edge_dt[is.na(edge_dt)] = 0
#     message('add_edge_width: Some NAs occur, unclear reason. For now NA->0.')
#
#     # Correct order is accomplished by the merge.
#     num_interacts = edge_dt[, num_interacts]
#     igraph::E(G)$width = rescale(sqrt(num_interacts), WIDTH_MIN, WIDTH_MAX)
#
#     return(G)
# }
# setup_vertices <- function(G, config) {
#     igraph::V(G)$vertex_types = infer_vertex_types(igraph::V(G)$name)
#     G = add_vertex_size(G)
#     return(G)
# }
# setup_edges <- function(G, config, use_adjpval) {
#     igraph::E(G)$arrow.size = config$EDGE_STYLE$ARROW_SIZE
#     igraph::E(G)$edge.color = purrr::map_chr(igraph::E(G)$edge_type, map_edgetype_to_color)
#     G = add_edge_width(G)
#     return(G)
# }
# setup_layout <- function(G, config, disperse) {
#     layout = igraph::layout_as_bipartite(
#         G,
#         types = igraph::V(G)$vertex_types,
#         hgap = config$LAYOUT$HGAP,
#         vgap = config$LAYOUT$VGAP
#     )
#     layout = layout[, 2:1]  # horizontal to vertical
#     layout = sort_layout_vertices(layout, G, config, disperse)
#     G$layout = layout
#     return(G)
# }
# rescale = function(v, min_, max_) {
#     ( v - min(v)) / ( max(v) - min(v) ) * (max_ - min_) + min_
# }
#
# #### Functions for layout ####
# get_edges_from_vertex = function(v_name, G) {
#     return(igraph::incident(G, v_name, mode='out'))
# }
#
# get_edges_to_vertex = function(v_name, G) {
#     return(igraph::incident(G, v_name, mode='in'))
# }
#
# count_up_and_down_edges = function(edges, config) {
#     color_up = config$EDGE_COLORING$COLOR_UP
#     color_down = config$EDGE_COLORING$COLOR_DOWN
#
#     counts = table(edges$edge.color)
#     num_up = as.integer(counts[color_up])
#     num_down = as.integer(counts[color_down])
#
#     num_up = ifelse(is.na(num_up), 0, num_up)
#     num_down = ifelse(is.na(num_down), 0, num_down)
#
#     return(list(N_UP=num_up, N_DOWN=num_down))
# }
#
# vertex_sort_key_atomic = function(v_name, from, G, config) {
#
#     if(from) {
#         edges = get_edges_from_vertex(v_name, G)
#     }
#     else {
#         edges = get_edges_to_vertex(v_name, G)
#     }
#     counts = count_up_and_down_edges(edges, config)
#     sort_key = as.integer(counts$N_UP - counts$N_DOWN)
#     return(sort_key)
# }
#
# vertex_sort_keys = function(v_names, from, G, config) {
#     sort_keys = purrr::map_dbl(
#         v_names,
#         ~ vertex_sort_key_atomic(.x, from, G, config)
#     )
#     return(sort_keys)
# }
# sort_layout_vertices = function(layout, G, config, disperse) {
#
#     vertex_names = igraph::V(G)$name
#     num_v = length(vertex_names)
#     midpoint = num_v/2
#     ligand_vertices = vertex_names[1:midpoint]
#     receptor_vertices = vertex_names[(midpoint+1):num_v]
#
#     ligand_keys = vertex_sort_keys(ligand_vertices, from=TRUE, G, config)
#     receptor_keys = vertex_sort_keys(receptor_vertices, from=FALSE, G, config)
#
#     # TODO: refactor into clear code
#     # The reorders layout rows to achieve sorting
#     layout_copy = copy(layout)
#     new_layout = copy(layout)
#     new_layout[order(ligand_keys), ] = layout_copy[1:midpoint, ]
#     new_layout[midpoint + order(receptor_keys), ] = layout_copy[(midpoint+1):num_v, ]
#
#     if(disperse) {
#         ## Make as function
#         # Get hgap
#         vgap = config$LAYOUT$HGAP
#         scale_factor = 3
#
#         new_layout[1:midpoint, 2][ligand_keys==0] = (
#             new_layout[1:midpoint, 2][ligand_keys==0] + scale_factor*vgap
#         )
#         new_layout[1:midpoint, 2][ligand_keys>0] = (
#             new_layout[1:midpoint, 2][ligand_keys>0] + 2*scale_factor*vgap
#         )
#
#         new_layout[(midpoint+1):num_v, 2][receptor_keys==0] = (
#             new_layout[(midpoint+1):num_v, 2][receptor_keys==0]
#             + scale_factor*vgap
#         )
#         new_layout[(midpoint+1):num_v, 2][receptor_keys>0] = (
#             new_layout[(midpoint+1):num_v, 2][receptor_keys>0]
#             + 2*scale_factor*vgap
#         )
#     }
#     return(new_layout)
# }
#
# #### Interactive networks ####
# extract_from_object = function(object) {
#     tryCatch(
#         error = function(cnd) {
#             return(list(
#                 tissue = object@parameters$object_name,
#                 ora = object@ora_default,
#                 dt_filtered = object@cci_detected
#             )
#             )
#         },
#         stop("Can't extract necessary data from object. Maybe check scDiffCom object.")
#     )
# }
# check_network_type_arg = function(network_type) {
#     if( !(network_type %in% c('bipartite', 'cells')) ) {
#         message = sprintf('network_type "%s" not in "bipartite", "cells"', network_type)
#         stop(message)
#     }
# }
# remove_LR_label = function(s) {
#     # Removes (L) or (R) label from celltype name.
#     # s: string vector
#     return(
#         sub(' [(][LR][)]', '', s)
#     )
# }
# map_bipartite_to_community = function(G) {
#
#     nodes = data.table::data.table(igraph::as_data_frame(G, what='vertices'))
#     edges = data.table::data.table(igraph::as_data_frame(G, what='edges'))
#
#     # Process edges
#     edges[, from := remove_LR_label(from)][,
#                                            to := remove_LR_label(to)]
#
#     # Process nodes
#     nodes[, name := remove_LR_label(name)]
#     nodes[, vertex_types := NULL]
#     nodes = nodes[, list(
#         vertex.size = mean(vertex.size),
#         vertex.counts = mean(vertex.counts)
#     ), by=name]
#
#     # Return new graph with similar properties for plotting
#     G_new = igraph::graph_from_data_frame(vertices=nodes, d=edges)
#     G_new$tissue = G$tissue
#     G_new$ora = G$ora
#     G_new$dt = G$dt
#     G_new$interacts_counts = G$interacts_counts
#     G_new$celltype_counts = G$celltype_counts
#
#     # New layout
#     # layout = igraph::layout_nicely(G_new, dim=2)  # determines best layout, likely calls fr
#     # layout = igraph::layout_with_fr(G_new)  # looks like default force-directed algorithm
#     layout = igraph::layout_with_sugiyama(G_new)$layout  # minimzes edge crossings
#     # layout = igraph::layout_in_circle(G_new)
#     G_new$layout = layout
#
#     return(G_new)
# }
# build_igraph = function(
#     ora,
#     dt_filtered,
#     tissue,
#     network_type
# ) {
#     ora_ct = ora$ER_CELLTYPES
#     ora_ct[, Tissue := tissue]
#     G = construct_graph(ora_ct, tissue, dt_filtered=dt_filtered)
#     G = setup_graph(G, use_adjpval=TRUE, disperse=TRUE)
#     if( network_type=='cells') {
#         G = map_bipartite_to_community(G)
#     }
#     return(G)
# }
# extract_node_edges = function(G, exclude_nonsign=TRUE) {
#     nodes = data.table::data.table(igraph::as_data_frame(G, what='vertices'))
#     edges = data.table::data.table(igraph::as_data_frame(G, what='edges'))
#
#     if( exclude_nonsign ){
#         # Will be moved into construct_graph
#         edges = edges[edge_type != 'none']
#     }
#
#     return(
#         list(nodes = nodes, edges = edges)
#     )
# }
# build_network_global = function(
#     nodes,
#     edges,
#     network_type
# ) {
#     nodes_ = data.table::copy(nodes)
#     edges_ = data.table::copy(edges)
#
#     title_template = '<p>%s</p><p>Median num cells:</p><p>Num interactions:</p>'
#     nodes_[,
#            'id' := name][,
#            'label' := name][,
#            'value' := vertex.size][,
#            'color' := 'green'][,  # Add: color.background, color.border, color.highlight.background/border, color.hover.background/hover
#            'shadow' := TRUE][,
#         # 'shapeProperties' := NULL][,
#            'title' := sprintf(title_template, name)]
#
#     if( network_type=='bipartite' ) {
#         nodes_[,
#                'group' := vertex_types][,
#                'level' := ifelse(vertex_types == TRUE, 1, 2)  # For hierarchical layout from vis
#                                         ]
#     } else if( network_type=='cells' ) {
#     } else {
#         stop('argument network_type not in "bipartite", "cells"')
#     }
#
#     edges_[,
#            color.color := edge.color][,
#            color.highlight := edge.color][,
#            color.hover := edge.color][,
#            smooth := TRUE]
#
#     return(
#         visNetwork::visNetwork(
#             nodes_, edges_,
#             width='100%',# height='100%',
#             main='Tissue cell:cell communication',
#             submain='...',
#             footer='Graph',
#             background='white'
#         )
#     )
# }
# get_network_components = function(
#     nodes,
#     edges,
#     layout,
#     network_type,
#     configure=FALSE
# ) {
#
#     if (configure) {
#         configure_component = . %>% visNetwork::visConfigure(
#             enabled=TRUE,
#             showButton=TRUE
#         )
#     } else {
#         configure_component = NULL
#     }
#
#     network_components = list(
#         network_global = build_network_global(
#             nodes, edges, network_type=network_type
#         ),
#         nodes_global = . %>% visNetwork::visNodes(
#             shape='dot',
#             physics=FALSE
#         ),
#         edges_global = . %>% visNetwork::visEdges(
#             shadow=TRUE,
#             arrows='middle',
#             smooth=list(enabled=TRUE, roundness=0.5)
#         ),
#         layout = . %>% visNetwork::visIgraphLayout(
#             layout='layout.norm',
#             layoutMatrix=layout
#         ),
#         options = . %>% visNetwork::visOptions(
#             highlightNearest=list(enabled=TRUE, degree=1, hover=TRUE),
#             selectedBy=list(variable="id", multiple=TRUE),
#             manipulation=TRUE
#         ),
#         interactive = . %>% visNetwork::visInteraction(
#             keyboard=list(enabled=TRUE),
#             multiselect=TRUE,
#             navigationButtons=TRUE
#         ),
#         configure = configure_component,
#         legend = . %>% visNetwork::visLegend(
#             enabled=TRUE,
#             useGroups=FALSE
#         )
#     )
#     return(network_components)
# }
# build_network = function(
#     network_global,
#     nodes_global=NULL,
#     edges_global=NULL,
#     layout=NULL,
#     legend=NULL,
#     options=NULL,
#     interactive=NULL,
#     configure=NULL
# ) {
#     vis_funcs = list(
#         nodes_global=nodes_global,
#         edges_global=edges_global,
#         layout=layout,
#         legend=legend,
#         options=options,
#         interactive=interactive,
#         configure=configure)
#     # For NULL arguments, use the identity as pipeline step
#     vis_funcs = purrr::map_if(vis_funcs, is.null, ~ . %>% identity())
#     return(
#         network_global %>%
#             (vis_funcs$nodes_global) %>%
#             (vis_funcs$edges_global) %>%
#             (vis_funcs$layout) %>%
#             (vis_funcs$legend) %>%
#             (vis_funcs$options) %>%
#             (vis_funcs$interactive) %>%
#             (vis_funcs$configure)
#     )
# }
# interactive_from_igraph = function(
#     G,
#     network_type,
#     exclude_nonsign=TRUE
# ) {
#
#     ne = extract_node_edges(G)
#     nodes = ne$nodes
#     edges = ne$edges
#
#     network_components = get_network_components(nodes, edges, G$layout, network_type)
#     interactive_network = do.call(build_network, network_components)
#     return(interactive_network)
# }
# build_interactive_network = function(
#     object,
#     network_type
# ){
#
#     check_network_type_arg(network_type)
#
#     extracted = extract_from_object(object)
#     tissue = extracted$tissue
#     ora = extracted$ora
#     dt_filtered = extracted$dt_filtered
#
#     G = build_igraph(ora, dt_filtered, tissue, network_type)
#     interactive_net = interactive_from_igraph(G, network_type)
#     return(interactive_net)
# }
#
# setup_graph_config = function() {
#     GRAPH_CONFIG = list(
#         EDGE_COLORING = list(
#             COLOR_UP = "#4285F4",  # blue
#             CUTOFF_UP = 1.1,
#             COLOR_DOWN = "#DB4437", # red
#             CUTOFF_DOWN = 1.1,
#             COLOR_BOTH = "#F4B400",
#             COLOR_ROBUST = "#0F9D58",
#             # COLOR_NONE = "#dde2e4"
#             COLOR_NONE = rgb(0.2, 0.2, 0.2, alpha=0.1)
#         ),
#
#         EDGE_STYLE = list(
#             ARROW_SIZE = 0.5,
#             WIDTH = 2.5  # will be adjusted by |OR|
#         ),
#
#         LAYOUT = list(
#             HGAP = 20,
#             VGAP = 20
#         ),
#
#         VERTEX_STYLE = list(
#             SIZE = 10,
#             LABEL_DIST = 1.5,
#             LABEL_CEX = 1.2,
#             COLOR = '#33FF66'
#         ),
#
#         LEGEND = list(
#             LEGEND_LABELS = c(
#                 'Significant, but small effect',
#                 'Upregulated',
#                 'Downregulated',
#                 'Altered'),
#             PCH = c(15),
#             CEX = 0.7,
#             PT.CEX = 1,
#             BG = '#CCCCCC',
#             NCOL = 2
#         )
#     )
#     return(GRAPH_CONFIG)
# }
