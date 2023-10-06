data(toy_data, package = 'CCPlotR')

test_that("Only accepts dataframe as cc_df", {
    expect_error(cc_network(seq_len(4)))
})

test_that("Checks for valid options", {
    expect_error(cc_network(toy_data, option = 1))
})

