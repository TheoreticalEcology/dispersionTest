library(DHARMa)

testData <- createData(overdispersion = 0,
                       sampleSize = 1000,
                       intercept = 0,
                       numGroups = 10,
                       randomEffectVariance = 0,
                       family = poisson())