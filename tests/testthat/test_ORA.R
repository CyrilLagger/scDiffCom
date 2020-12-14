context("ORA")

object = readRDS("../../../data_scAgeCom/testing/scDiffCom_test_facs_liver.rds")

test_that("Component functions return", {
    ora_go = scDiffCom::build_ora_dt(
        object@cci_detected,
        object@parameters$threshold_logfc,
        'UP',
        'GO_TERMS',
        'mouse'
    )
    expect_equal(dim(ora_go), c(2687, 13))
})