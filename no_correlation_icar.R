require(data.table)  # for data mangling
require(tidyverse)  # for data mangling
require(ggplot2)  # for plotting
require(hexbin)  # for plotting
require(bayesplot)  # for plotting Stan outputs
require(knitr)  # for Rmarkdown
require(kableExtra)  # for Rmarkdown
require(cmdstanr)  # for Stan
library(dplyr)
library(stringr)
library(glue)
library(posterior)

# set colour scheme for Stan
bayesplot::color_scheme_set("brewer-RdYlBu")

# change this line as needed to point to the directory where all output from
# this lab is stored
out.dir <- "/Users/44775/Documents"


ws10 <- read.csv("pre_icar_data/y_ws10.csv", header=FALSE)
tp <- read.csv("pre_icar_data/y_tp.csv", header = FALSE)

node1 <- read.csv("pre_icar_data/node1.csv", header = FALSE)
node2 <- read.csv("pre_icar_data/node2.csv", header = FALSE)

node1 <- as.integer(node1[[1]])  # Convert node1 to a vector
node2 <- as.integer(node2[[1]])  # Convert node2 to a vector
y1 <- as.vector(ws10[[1]])  # Convert ws10 to a numeric vector
y2 <- as.vector(tp[[1]])  # Convert ws10 to a numeric vector

# Data list
data = list(
  N = length(y1), 
  N_edges = length(node1), 
  node1 = node1, 
  node2 = node2,
  y1 = (y1 - mean(y1))/sd(y1),
  y2 = (y2 - mean(y2))/sd(y2)
)

model1_text <- "
data {
  int<lower=0> N;
  int<lower=0> N_edges;
  array[N_edges] int<lower=1, upper=N> node1; // node1[i] adjacent to node2[i]
  array[N_edges] int<lower=1, upper=N> node2; // and node1[i] < node2[i]
  
  array[N] real y1; // var1 outcome
  array[N] real y2; // var2 outcome
}
parameters {
  real beta0_1; // intercept
  real beta0_2;
  real<lower=0.00001> sigma_1; // overall standard deviation
  real<lower=0.00001> sigma_2;
  
  vector[2*N] phi; // spatial effects
  real<lower=0.00001> sigma_phi1;
  real<lower=0.00001> sigma_phi2;
}
transformed parameters {
  vector[N] phi1 = phi[1:N];
  vector[N] phi2 = phi[(N+1):(2*N)];
}
model {
  // Model with phi1 and phi2
  y1 ~ normal(beta0_1 + phi1 * sigma_phi1, sigma_1);
  y2 ~ normal(beta0_2 + phi2 * sigma_phi2, sigma_2);

  beta0_1 ~ normal(0, 5);
  beta0_2 ~ normal(0, 5);
  
  sigma_1 ~ normal(0.6, 0.03);  // Ensure sigma_1 is not zero
  sigma_2 ~ normal(0.6, 0.03);
  
  sigma_phi1 ~ normal(0.4, 0.02);
  sigma_phi2 ~ normal(0.4, 0.02);
  
  target += -0.5 * (dot_self(phi1[node1] - phi1[node2]) + dot_self(phi2[node1] - phi2[node2]));
  
  sum(phi1) ~ normal(0, 0.001 * N);
  sum(phi2) ~ normal(0, 0.001 * N);
}

generated quantities {
  vector[N] mu1 = beta0_1 + phi1 * sigma_phi1;
  vector[N] mu2 = beta0_2 + phi2 * sigma_phi2;
}
"

# cmdstanr requires the model to be written to a file
model1_filename <- cmdstanr::write_stan_file(
  gsub('\t',' ',model1_text),
  dir = out.dir,
  basename = NULL,
  force_overwrite = FALSE,
  hash_salt = ""
)

# compile Stan model
model1_compiled_cmdstanr <- cmdstanr::cmdstan_model(model1_filename)

# Fit the model
model1_fit <- model1_compiled_cmdstanr$sample(
  data = data,
  seed = 123,
  chains = 2,
  parallel_chains = 2,
  refresh = 500,
  iter_warmup = 1000,
  iter_sampling = 3000,
  save_warmup = TRUE
)

# save output to RDS
model1_fit$save_object(file = file.path(out.dir, "simple_ICAR_compiled_cmdstanr.rds"))

# load output from RDS
model1_fit <- readRDS(file.path(out.dir, "simple_ICAR_compiled_cmdstanr.rds"))
model1_pars <- c("beta0_1", "beta0_2", "sigma_1", "sigma_2", "phi1[1]", "mu1[1]", "phi2[1]", "mu2[1]", "mu2[2]", "sigma_phi1", "sigma_phi2")


model1_fit$summary(
  variables = model1_pars,
  posterior::default_summary_measures(),
  posterior::default_convergence_measures(),
  extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975)))


# Save output
mu1 <- c()
mu2 <- c()
N <- length(y1)

# Loop over each node index (based on the length of node1)
for (i in 1:N) {
  
  # Get the summary statistics for phi[i] (adjust indexing to R style)
  stat1 <- model1_fit$summary(variables = paste("mu1[", i, "]", sep = ""))
  stat2 <- model1_fit$summary(variables = paste("mu2[", i, "]", sep = ""))
  
  # Add the mean of phi[i] to the sum
  mu1 <- append(mu1, stat1$median*sd(y1) + mean(y1))
  mu2 <- append(mu2, stat2$median*sd(y2) + mean(y2))


write.csv(data.frame(mu1), "post_icar_data/one_variable/mu1.csv", row.names = FALSE)

write.csv(data.frame(mu2), "post_icar_data/one_variable/mu2.csv", row.names = FALSE)