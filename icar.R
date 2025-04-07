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


#Data
ws10 <- read.csv("pre_icar_data/y_ws10.csv", header=FALSE)
tp <- read.csv("pre_icar_data/y_tp.csv", header = FALSE)
corr <- read.csv("pre_icar_data/corr.csv", header = FALSE)

node1 <- read.csv("pre_icar_data/node1.csv", header = FALSE)
node2 <- read.csv("pre_icar_data/node2.csv", header = FALSE)

rho_node1 <- read.csv("pre_icar_data/rho_node1.csv", header = FALSE)
rho_node2 <- read.csv("pre_icar_data/rho_node2.csv", header = FALSE)

node1 <- as.integer(node1[[1]])  # Convert node1 to a vector
node2 <- as.integer(node2[[1]])  # Convert node2 to a vector

y1 <- as.vector(ws10[[1]])  # Convert ws10 to a numeric vector
y2 <- as.vector(tp[[1]])  # Convert ws10 to a numeric vector

theta_corr <- as.vector(corr[[1]])


# Data list
data = list(
  N = length(y1), 
  N_edges = length(node1),
  node1 = node1, 
  node2 = node2,
  y1 = (y1-mean(y1))/sd(y1),
  y2 = (y2-mean(y2))/sd(y2),
  theta_corr = theta_corr
)


model1_text <-"
data {

  int<lower=0> N;              
  int<lower=0> N_edges;         
  array[N_edges] int node1;     
  array[N_edges] int node2;     

  vector[N] y1;        
  vector[N] y2;
  
  vector[N] theta_corr;
}

parameters {

  real beta0_1;
  real beta0_2;
  real<lower=0> sigma_1;
  real<lower=0> sigma_2;

  vector[2*N] phi; 
  real<lower=0> sigma_phi1;
  real<lower=0> sigma_phi2;
  
  vector<lower=-1, upper=1>[N] theta;
}

transformed parameters {
  vector[N] phi1 = phi[1:N];
  vector[N] phi2 = phi[(N+1):(2*N)];
  
  vector[N_edges] theta_edge = (theta[node1] + theta[node2])/2;
}

model {

  y1 ~ normal(beta0_1 + phi1 * sigma_phi1, sigma_1 + 1e-6);
  y2 ~ normal(beta0_2 + phi2 * sigma_phi2, sigma_2 + 1e-6);

  beta0_1 ~ normal(0,5);
  beta0_2 ~ normal(0,5);
  
  sigma_1 ~ normal(0.6, 0.03);  // Ensure sigma_1 is not zero
  sigma_2 ~ normal(0.6, 0.03);
  
  sigma_phi1 ~ normal(0.4, 0.02);
  sigma_phi2 ~ normal(0.4, 0.02);
  
  theta ~ normal(theta_corr, 0.5);
  
  target += -0.5 * (dot_self(phi1[node1] - phi1[node2]) + dot_self(phi2[node1] - phi2[node2]));
  
  target += dot_product( phi1 .* phi2, theta);
  target += dot_product( phi1[node1] .* phi2[node2] + phi2[node1] .* phi1[node2], theta_edge);
  
  target += -1 * dot_product(phi1 .* phi1 + phi2 .* phi2, abs(theta));
  target += -1 * dot_product(phi1[node1] .* phi1[node1] + phi2[node1] .* phi2[node1], abs(theta_edge));
  
  target += -0.5 * dot_self(theta[node1] - theta[node2]);

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


# Fit the model with the updated init
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
model1_pars <- c("lp__","beta0_1", "beta0_2", "sigma_1", "sigma_2", "mu1[1]", "mu1[2]", "mu2[1]", "mu2[2]", "sigma_phi1","sigma_phi2", "theta[1]")


model1_fit$summary(
  variables = model1_pars,
  posterior::default_summary_measures(),
  posterior::default_convergence_measures(),
  extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975)))


# Trace Plots to check convergence

all_pars <- c("lp__","beta0_1", "beta0_2", "sigma_1", "sigma_2", "phi[1]", "sigma_phi1","sigma_phi2", "theta[1]")

all_draws <- model1_fit$draws(
  variables = c("lp__",all_pars),
  inc_warmup = FALSE,
  format = "draws_array"
)
all_draws_df <- as_draws_df(all_draws)
write_csv(all_draws_df, glue("{out.dir}/R/all_draws.csv"))

### Test parameters
p <- bayesplot:::mcmc_trace(all_draws,  
                            pars = all_pars, 
                            n_warmup = 0,
                            facet_args = list(nrow = 5)
)
# use this for png
ggsave(file = glue("{out.dir}/R/all_pars_traceplot.png"), 
       plot = p, 
       h = 12, 
       w = 8
)


# Save results
mu1 <- c()
mu2 <- c()
N <- length(y1)

# Loop over each node index (based on the length of node1)
for (i in 1:N) {
  
  # Get the summary statistics for phi[i] (adjust indexing to R style)
  stat1 <- model1_fit$summary(variables = paste("mu1[", i, "]", sep = ""))
  stat2 <- model1_fit$summary(variables = paste("mu2[", i, "]", sep = ""))
  
  # Add the mean of phi[i] to the sum
  mu1 <- append(mu1, stat1$median * sd(y1) + mean(y1))
  mu2 <- append(mu2, stat2$median * sd(y2) + mean(y2))
}

correlations <- c()

# Loop over each node index (based on the length of node1)
for (i in 1:N) {
  
  # Get the summary statistics for phi[i] (adjust indexing to R style)
  stat1 <- model1_fit$summary(variables = paste("theta[", i, "]", sep = ""))
  
  # Add the mean of phi[i] to the sum
  correlations <- append(correlations, stat1$median)
}


write.csv(data.frame(correlations), "post_icar_data/two_variable/correlations.csv", row.names = FALSE)

write.csv(data.frame(mu1), "post_icar_data/two_variable/mu1.csv", row.names = FALSE)

write.csv(data.frame(mu2), "post_icar_data/two_variable/mu2.csv", row.names = FALSE)