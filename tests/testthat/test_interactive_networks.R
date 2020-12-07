context("Interactive networks")

object = readRDS("../../../data_scAgeCom/testing/scDiffCom_test_facs_liver.rds")
G = construct_graph(object)

test_that("process_celltype_pairs_enrichment", {
    res = process_celltype_pairs_enrichment(object@ora_default$ER_CELLTYPES)
    expect_equal(dim(res), c(49, 16))
    expect_true( !('tissue' %in% tolower(names(res))) )
    expect_equal(
        sum(object@cci_detected$ER_CELLTYPES == 'natural killer cell_natural killer cell'),
        res[Ligand_cell == 'natural killer cell' & Receptor_cell == 'natural killer cell', COUNTS_TOTAL]
    )
})

test_that("average_celltype_counts", {
    counts = average_celltype_counts(object)
    expect_equal(dim(counts), c(7, 2))
    expect_equal(sum(counts[, 2]), 1429.5)
})

test_that("num_interactions_cellpairs", {
    counts = num_interactions_cellpairs(object)
    expect_equal(dim(counts), c(49, 3))
    expect_equal(sum(counts[, 3]), 4003)
})

test_that("construct_graph", {
    expect_equal(length(G), 10)
    expect_equal(length(igraph::V(G)), 14)
    expect_equal(length(igraph::E(G)), 49)
})

test_that("add_vertex_size", {
    g = add_vertex_size(G)
    expect_equal(round(sum(igraph::V(g)$vertex.size)), 139)
    expect_equal(max(igraph::V(g)$vertex.size), 20)
})

test_that("add_edge_width", {
    g = add_edge_width(G)
    expect_equal(max(igraph::E(g)$width), 10)
    expect_equal(min(igraph::E(g)$width), 0)
})