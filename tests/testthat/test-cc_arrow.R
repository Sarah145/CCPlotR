data(toy_data, package = 'CCPlotR')

test_that("Only accepts dataframe as cc_df", {
    expect_error(cc_arrow(seq_len(4)))
})

test_that("Checks for valid options", {
    expect_error(cc_arrow(toy_data, option = 'C'))
})

test_that("Checks cell_types vector", {
    expect_error(cc_arrow(toy_data, cell_types = seq_len(4)))
})
