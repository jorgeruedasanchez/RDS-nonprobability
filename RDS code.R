# Codes for generation of networks, RDS estimators, variance estimation  and two types of #Confidence Intervals

# Generation of networks

library(igraph)

population_size <- 10000
small_world_nei <- 5

# Barabási and Albert (1999) Scale-free network
network1 <- barabasi.game(n = population_size, directed = FALSE)

# Watts and Strogatz (1998) Small-world network low p 
network2 <- watts.strogatz.game(dim = 1, size = population_size, nei = small_world_nei, p = 0.01)

# Watts and Strogatz (1998) Small-world network high p
network3 <- watts.strogatz.game(dim = 1, size = population_size, nei = small_world_nei, p = 0.5)


# Computation of network metrics

avg_deg <- mean(degree(network1))
max_deg <- max(degree(network1))
avg_path <- average.path.length(network1)
clustering <- transitivity(network1, type = "global")

metrics <- data.frame(
  AverageDegree = round(avg_deg, 2),
  MaxDegree = max_deg,
  AvgPathLength = round(avg_path, 2),
  ClusteringCoeff = round(clustering, 2)
)

print(metrics)


# Computation of RDS estimators

library(RDS)

rds_data <- as.rds.data.frame(dat, id = "id", recruiter.id = "recruiter.id", network.size = "degree")

# RDS-I estimator 
rdsi <- RDS.I.estimates(rds_data, outcome.variable = "y")
esti <- rdsi$estimate

print(esti)

# RDS-II estimator 
rdsii <- RDS.II.estimates(rds_data, outcome.variable = "y")
estii <- rdsii$estimate

print(estii)



# Gile's Successive Sampling RDS-SS
rds_ss <- RDS.SS.estimates(rds_data, outcome.variable = "y", N = 10000)
estss <- rds_ss$estimate

print(estss)


# Variance estimation and Confidence Intervals for RDS estimators

library(RDS)

# RDS data frame conversion
rds_data <- as.rds.data.frame(dat, id = "id", recruiter.id = "recruiter.id", network.size = "degree")


# Salganik Bootstrap (Sal-BS) (Salganik, 2006) Variance estimation and Quantile Confidence  #Interval for RDS-I estimator
CI1 <- RDS.bootstrap.intervals(
  rds.data = rds_data, 
  weight.type = "RDS-I",
  confidence.level = 0.95,
  outcome.variable = "y",
  ci.type="percentile", 
  number.of.bootstrap.samples = 100
)

print(CI1)

# Salganik Bootstrap (Sal-BS) (Salganik, 2006) Variance estimation and Quantile Confidence #Interval for RDS-II estimator
CI2 <- RDS.bootstrap.intervals(
  rds.data = rds_data,
  weight.type = "RDS-II",
  confidence.level = 0.95,
  outcome.variable = "y",
  ci.type="percentile", 
  number.of.bootstrap.samples = 100
)

print(CI2)

# Successive Sampling Bootstrap (SS-BS) Variance estimation and Studentized Confidence Interval #for Successive Sampling (RDS-SS) estimator
CI3 <- RDS.bootstrap.intervals(
  rds.data = rds_data,
  weight.type = "Gile's SS", 
  confidence.level = 0.95,    
  outcome.variable = "y",
  N = population_size,
  ci.type="t",
  number.of.bootstrap.samples = 100
)

print(CI3)

