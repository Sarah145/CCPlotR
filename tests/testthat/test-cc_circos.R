data(toy_data, package = 'CCPlotR')

test_that("Only accepts dataframe as cc_df", {
    expect_error(cc_circos(seq_len(4)))
})

test_that("Checks for valid options", {
    expect_error(cc_circos(toy_data, option = 1))
})

test_that("Only accepts dataframe as exp_df", {
    expect_error(cc_circos(toy_data, option = 'C', exp_df = seq_len(4)))
})

