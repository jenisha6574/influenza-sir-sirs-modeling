library(tidyverse)
library(odin)

flu_model_vax <- odin({
  deriv(S) <- -beta * S * I + omega * R - v * S
  deriv(I) <-  beta * S * I - gamma * I
  deriv(R) <-  gamma * I - omega * R + v * S
  
  initial(S) <- 9990
  initial(I) <-   10
  initial(R) <-    0
  
  beta  <- user()
  gamma <- user()
  omega <- user()
  v     <- user()
})


gamma_val <- 0.2          
time      <- seq(0, 200, by = 1)

run_scenario <- function(beta, omega, v, label) {
  model <- flu_model_vax$new(
    beta  = beta,
    gamma = gamma_val,
    omega = omega,
    v     = v
  )
  
  as_tibble(as.data.frame(model$run(time))) %>%
    mutate(scenario = label)
}

beta_low  <- 0.0001
beta_high <- 0.0005
omega_SIR <- 0
v_none    <- 0

sir_low_beta <- run_scenario(beta_low,  omega_SIR, v_none,
                             paste("SIR: low beta =", beta_low))
sir_high_beta <- run_scenario(beta_high, omega_SIR, v_none,
                              paste("SIR: high beta =", beta_high))

sir_beta_compare <- bind_rows(sir_low_beta, sir_high_beta) %>%
  select(t, I, scenario)

fig_SIR_beta <- sir_beta_compare %>%
  ggplot(aes(x = t, y = I, color = scenario)) +
  geom_line(linewidth = 1) +
  labs(
    title = paste(
      "SIR Model: Effect of beta | gamma =", gamma_val,
      "| omega = 0 | vaccine rate = 0"
    ),
    x = "Time (days)",
    y = "Number of infected people",
    color = "Scenario"
  ) +
  theme_classic()

fig_SIR_beta

omega_SIRS <- 0.05   # immunity loss = every 20 days

sirs_low_beta <- run_scenario(beta_low,  omega_SIRS, v_none,
                              paste("SIRS: low beta =", beta_low))
sirs_high_beta <- run_scenario(beta_high, omega_SIRS, v_none,
                               paste("SIRS: high beta =", beta_high))

sirs_beta_compare <- bind_rows(sirs_low_beta, sirs_high_beta) %>%
  select(t, I, scenario)

fig_SIRS_beta <- sirs_beta_compare %>%
  ggplot(aes(x = t, y = I, color = scenario)) +
  geom_line(linewidth = 1) +
  labs(
    title = paste(
      "SIRS Model: Effect of beta | gamma =", gamma_val,
      "| omega =", omega_SIRS, "| vaccine rate = 0"
    ),
    x = "Time (days)",
    y = "Number of infected people",
    color = "Scenario"
  ) +
  theme_classic()

fig_SIRS_beta

##  SIRS, low vs high omega 
beta_fixed <- 0.0003
omega_low  <- 0.02    # this means immunity lasts longer
omega_high <- 0.08    # this meansimmunity fades faster

sirs_low_omega <- run_scenario(beta_fixed, omega_low, v_none,
                               paste("SIRS: low omega =", omega_low))
sirs_high_omega <- run_scenario(beta_fixed, omega_high, v_none,
                                paste("SIRS: high omega =", omega_high))

sirs_omega_compare <- bind_rows(sirs_low_omega, sirs_high_omega) %>%
  select(t, I, scenario)

fig_SIRS_omega <- sirs_omega_compare %>%
  ggplot(aes(x = t, y = I, color = scenario)) +
  geom_line(linewidth = 1) +
  labs(
    title = paste(
      "SIRS Model: Effect of immunity loss omega",
      "| beta =", beta_fixed,
      "| gamma =", gamma_val,
      "| vaccine rate = 0"
    ),
    x = "Time (days)",
    y = "Number of infected people",
    color = "Scenario"
  ) +
  theme_classic()

fig_SIRS_omega
beta_vax   <- 0.0003
omega_SIR  <- 0
omega_SIRS <- 0.05
v_low      <- 0.01
v_high     <- 0.5      # strong campaign

##SIR, no vaccine vs vaccine 
sir_no_vax <- run_scenario(beta_vax, omega_SIR, v_low,
                           "SIR: no vaccine")
sir_with_vax <- run_scenario(beta_vax, omega_SIR, v_high,
                             paste("SIR: vaccine rate =", v_high))

sir_vax_compare <- bind_rows(sir_no_vax, sir_with_vax) %>%
  select(t, I, scenario)

fig_SIR_vax <- sir_vax_compare %>%
  ggplot(aes(x = t, y = I, color = scenario)) +
  geom_line(linewidth = 1) +
  labs(
    title = paste(
      "SIR Model: Effect of vaccination",
      "| beta =", beta_vax,
      "| gamma =", gamma_val,
      "| omega = 0"
    ),
    x = "Time (days)",
    y = "Number of infected people",
    color = "Scenario"
  ) +
  theme_classic()

fig_SIR_vax

##  no vaccine vs vaccine 
sirs_no_vax <- run_scenario(beta_vax, omega_SIRS, v_low,
                            "SIRS: no vaccine")
sirs_with_vax <- run_scenario(beta_vax, omega_SIRS, v_high,
                              paste("SIRS: vaccine rate =", v_high))

sirs_vax_compare <- bind_rows(sirs_no_vax, sirs_with_vax) %>%
  select(t, I, scenario)

fig_SIRS_vax <- sirs_vax_compare %>%
  ggplot(aes(x = t, y = I, color = scenario)) +
  geom_line(linewidth = 1) +
  labs(
    title = paste(
      "SIRS Model: Effect of vaccination",
      "| beta =", beta_vax,
      "| gamma =", gamma_val,
      "| omega =", omega_SIRS
    ),
    x = "Time (days)",
    y = "Number of infected people",
    color = "Scenario"
  ) +
  theme_classic()

fig_SIRS_vax
