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
sp <- read.csv("pre_icar_data/y_sp.csv", header = FALSE)

rho_corr <- read.csv("pre_icar_data/corr.csv", header = FALSE)
gam_corr <- read.csv("pre_icar_data/gam_corr.csv", header = FALSE)
del_corr <- read.csv("pre_icar_data/del_corr.csv", header = FALSE)

node1 <- read.csv("pre_icar_data/node1.csv", header = FALSE)
node2 <- read.csv("pre_icar_data/node2.csv", header = FALSE)

node1 <- as.integer(node1[[1]])  # Convert node1 to a vector
node2 <- as.integer(node2[[1]])  # Convert node2 to a vector

y1 <- as.vector(ws10[[1]])  # Convert ws10 to a numeric vector
y2 <- as.vector(tp[[1]])  # Convert ws10 to a numeric vector
y3 <- as.vector(sp[[1]])

rho_corr <- as.vector(rho_corr[[1]])
gam_corr <- as.vector(gam_corr[[1]])
del_corr <- as.vector(del_corr[[1]])


# Data list
data = list(
  N = length(y1), 
  N_edges = length(node1),
  node1 = node1, 
  node2 = node2,
  y1 = (y1-mean(y1))/sd(y1),
  y2 = (y2-mean(y2))/sd(y2),
  y3 = (y3-mean(y3))/sd(y3),
  rho_corr = rho_corr,
  gam_corr = gam_corr,
  del_corr = del_corr
)


model1_text <-"
functions {
  real interaction_term(vector phi1, vector phi2, array[] int node1, array[] int node2, vector theta, vector theta_edge) {
    real out = 0;
    
    out += dot_product(phi1 .* phi2, theta);
    out += dot_product(phi1[node1] .* phi2[node2] + phi2[node1] .* phi1[node2], theta_edge);

    out += -1 * dot_product(phi1 .* phi1 + phi2 .* phi2, abs(theta));
    out += -1 * dot_product(phi1[node1] .* phi1[node1] + phi2[node1] .* phi2[node1], abs(theta_edge));
    
    return out;
  }
}


data {

  int<lower=0> N;              
  int<lower=0> N_edges;         
  array[N_edges] int node1;     
  array[N_edges] int node2;     

  vector[N] y1;        
  vector[N] y2;
  vector[N] y3;
  
  vector[N] rho_corr;
  vector[N] gam_corr;
  vector[N] del_corr;
}

parameters {

  real beta0_1;
  real beta0_2;
  real beta0_3;
  
  real<lower=0> sigma_1;
  real<lower=0> sigma_2;
  real<lower=0> sigma_3;

  vector[3*N] phi; 
  real<lower=0> sigma_phi1;
  real<lower=0> sigma_phi2;
  real<lower=0> sigma_phi3;
  
  vector<lower=-1, upper=1>[N] rho;
  vector<lower=-1, upper=1>[N] gamma;
  vector<lower=-1, upper=1>[N] delta;
}

transformed parameters {
  vector[N] phi1 = phi[1:N];
  vector[N] phi2 = phi[(N+1):(2*N)];
  vector[N] phi3 = phi[(2*N+1):(3*N)];
  
  vector[N_edges] rho_edge = (rho[node1] + rho[node2])/2;
  vector[N_edges] gamma_edge = (gamma[node1] + gamma[node2])/2;
  vector[N_edges] delta_edge = (delta[node1] + delta[node2])/2;
}

model {

  y1 ~ normal(beta0_1 + phi1 * sigma_phi1, sigma_1 + 1e-6);
  y2 ~ normal(beta0_2 + phi2 * sigma_phi2, sigma_2 + 1e-6);
  y3 ~ normal(beta0_3 + phi3 * sigma_phi3, sigma_3 + 1e-6);

  beta0_1 ~ normal(0,5);
  beta0_2 ~ normal(0,5);
  beta0_3 ~ normal(0,5);
  
  sigma_1 ~ normal(0.6, 0.03);
  sigma_2 ~ normal(0.6, 0.03);
  sigma_3 ~ normal(0.6, 0.03);
  
  sigma_phi1 ~ normal(0.4, 0.02);
  sigma_phi2 ~ normal(0.4, 0.02);
  sigma_phi3 ~ normal(0.4, 0.02);
  
  rho ~ normal(0.5, 1);
  gamma ~ normal(0, 1);
  delta ~ normal(-0.5, 1);
  
  target += -0.5 * (dot_self(phi1[node1] - phi1[node2]) + dot_self(phi2[node1] - phi2[node2]) + dot_self(phi3[node1] - phi3[node2]));
  
  target += interaction_term(phi1, phi2, node1, node2, rho, rho_edge);
  target += interaction_term(phi1, phi3, node1, node2, gamma, gamma_edge);
  target += interaction_term(phi2, phi3, node1, node2, delta, delta_edge);
  
  target += -0.5 * (dot_self(rho[node1] - rho[node2]) + dot_self(gamma[node1] - gamma[node2]) + dot_self(delta[node1] - delta[node2]));

  sum(phi1) ~ normal(0, 0.001 * N);
  sum(phi2) ~ normal(0, 0.001 * N);
  sum(phi3) ~ normal(0, 0.001 * N);
  
  sd(rho) ~ normal(1,0.015);
  sd(gamma) ~ normal(1,0.015);
  sd(delta) ~ normal(1,0.015);
  
}

generated quantities {

  vector[N] mu1 = beta0_1 + phi1 * sigma_phi1;
  vector[N] mu2 = beta0_2 + phi2 * sigma_phi2;
  vector[N] mu3 = beta0_3 + phi3 * sigma_phi3;
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



mu1 <- c()
mu2 <- c()
mu3 <- c()
N <- length(y1)

# Loop over each node index (based on the length of node1)
for (i in 1:N) {
  
  # Get the summary statistics for phi[i] (adjust indexing to R style)
  stat1 <- model1_fit$summary(variables = paste("mu1[", i, "]", sep = ""))
  stat2 <- model1_fit$summary(variables = paste("mu2[", i, "]", sep = ""))
  stat3 <- model1_fit$summary(variables = paste("mu3[", i, "]", sep = ""))
  
  # Add the mean of phi[i] to the sum
  mu1 <- append(mu1, stat1$median * sd(y1) + mean(y1))
  mu2 <- append(mu2, stat2$median * sd(y2) + mean(y2))
  mu3 <- append(mu3, stat3$median * sd(y3) + mean(y3))
}


correlations_rho <- c()
correlations_gam <- c()
correlations_del <- c()
N <- length(y1)

# Loop over each node index (based on the length of node1)
for (i in 1:N) {
  
  # Get the summary statistics for phi[i] (adjust indexing to R style)
  stat1 <- model1_fit$summary(variables = paste("rho[", i, "]", sep = ""))
  stat2 <- model1_fit$summary(variables = paste("gamma[", i, "]", sep = ""))
  stat3 <- model1_fit$summary(variables = paste("delta[", i, "]", sep = ""))
  
  # Add the mean of phi[i] to the sum
  correlations_rho <- append(correlations_rho, stat1$median)
  correlations_gam <- append(correlations_gam, stat2$median)
  correlations_del <- append(correlations_del, stat3$median)
}


write.csv(data.frame(correlations_rho), "post_icar_data/three_variable/correlations_rho.csv", row.names = FALSE)
write.csv(data.frame(correlations_gam), "post_icar_data/three_variable/correlations_gam.csv", row.names = FALSE)
write.csv(data.frame(correlations_del), "post_icar_data/three_variable/correlations_del.csv", row.names = FALSE)

write.csv(data.frame(mu1), "post_icar_data/three_variable/mu1.csv", row.names = FALSE)
write.csv(data.frame(mu2), "post_icar_data/three_variable/mu2.csv", row.names = FALSE)
write.csv(data.frame(mu3), "post_icar_data/three_variable/mu3.csv", row.names = FALSE)