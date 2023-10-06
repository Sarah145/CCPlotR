data(toy_data, package = 'CCPlotR')

test_that("Only accepts dataframe as cc_df", {
    expect_error(cc_sigmoid(seq_len(4)))
})
