phy = ape::rtree(n = 300)

test_that("vcv should have the same results with ape::vcv", {
    expect_equal(ape::vcv(phy, corr = F), vcv2(phy, corr = F))
    expect_equal(ape::vcv(phy, corr = T), vcv2(phy, corr = T))
})

# microbenchmark::microbenchmark(ape::vcv(phy, corr = F), vcv2(phy, corr = F),
# times = 10) microbenchmark::microbenchmark(ape::vcv(phy, corr = T), vcv2(phy,
# corr = T), times = 10)
