# Définition des fonctions pour les dérivées du modèle SIR
dSIRdt <- function(S, I, beta, gamma, N) {
  dSdt = -beta * S * I / N
  dIdt = beta * S * I / N - gamma * I
  dRdt = gamma * I
  return(c(dSdt, dIdt, dRdt))
}

# Définition des fonctions pour les dérivées du modèle SIS
dSISdt <- function(S, I, beta, gamma, N) {
  dSdt = -beta * S * I / N + gamma * I
  dIdt = beta * S * I / N - gamma * I
  return(c(dSdt, dIdt))
}
# Fonction simplifiée pour estimer l'erreur entre deux pas de temps
estimate_error <- function(step1, step2) {
  max(abs(step1 - step2))
}

# Fonction RK6_step 
RK6_step <- function(model_type, S, I, beta, gamma, N, dt) {
  # Les calculs RK6 communs
  a21 = 1/3
  a31 = 1/12; a32 = 1/12
  a41 = 1/8; a42 = 0; a43 = 3/8
  a51 = 13/54; a52 = 0; a53 = -27/54; a54 = 42/54
  a61 = 1/6; a62 = 0; a63 = 2/3; a64 = -1/6; a65 = 1/6
  c1 = 1/6; c2 = 1/3; c3 = 1/3; c4 = 1/6; c5 = 0; c6 = 0
  
  # k1 calculations
  if (model_type == "SIR") {
    derivs = dSIRdt(S, I, beta, gamma, N)
  } else if (model_type == "SIS") {
    derivs = dSISdt(S, I, beta, gamma, N)
  }
  k1 = dt * derivs
  
  # k2 calculations
  if (model_type == "SIR") {
    derivs = dSIRdt(S + a21 * k1[1], I + a21 * k1[2], beta, gamma, N)
  } else if (model_type == "SIS") {
    derivs = dSISdt(S + a21 * k1[1], I + a21 * k1[2], beta, gamma, N)
  }
  k2 = dt * derivs
  
  # k3 calculations
  if (model_type == "SIR") {
    derivs = dSIRdt(S + a31 * k1[1] + a32 * k2[1], I + a31 * k1[2] + a32 * k2[2], beta, gamma, N)
  } else if (model_type == "SIS") {
    derivs = dSISdt(S + a31 * k1[1] + a32 * k2[1], I + a31 * k1[2] + a32 * k2[2], beta, gamma, N)
  }
  k3 = dt * derivs
  
  # k4 calculations
  if (model_type == "SIR") {
    derivs = dSIRdt(S + a41 * k1[1] + a42 * k2[1] + a43 * k3[1], I + a41 * k1[2] + a42 * k2[2] + a43 * k3[2], beta, gamma, N)
  } else if (model_type == "SIS") {
    derivs = dSISdt(S + a41 * k1[1] + a42 * k2[1] + a43 * k3[1], I + a41 * k1[2] + a42 * k2[2] + a43 * k3[2], beta, gamma, N)
  }
  k4 = dt * derivs
  
  # k5 calculations
  if (model_type == "SIR") {
    derivs = dSIRdt(S + a51 * k1[1] + a52 * k2[1] + a53 * k3[1] + a54 * k4[1], I + a51 * k1[2] + a52 * k2[2] + a53 * k3[2] + a54 * k4[2], beta, gamma, N)
  } else if (model_type == "SIS") {
    derivs = dSISdt(S + a51 * k1[1] + a52 * k2[1] + a53 * k3[1] + a54 * k4[1], I + a51 * k1[2] + a52 * k2[2] + a53 * k3[2] + a54 * k4[2], beta, gamma, N)
  }
  k5 = dt * derivs
  
  # k6 calculations
  if (model_type == "SIR") {
    derivs = dSIRdt(S + a61 * k1[1] + a62 * k2[1] + a63 * k3[1] + a64 * k4[1] + a65 * k5[1], I + a61 * k1[2] + a62 * k2[2] + a63 * k3[2] + a64 * k4[2] + a65 * k5[2], beta, gamma, N)
  } else if (model_type == "SIS") {
    derivs = dSISdt(S + a61 * k1[1] + a62 * k2[1] + a63 * k3[1] + a64 * k4[1] + a65 * k5[1], I + a61 * k1[2] + a62 * k2[2] + a63 * k3[2] + a64 * k4[2] + a65 * k5[2], beta, gamma, N)
  }
  k6 = dt * derivs
  
  # Combine k-values for final result
  S_next = S + c1 * k1[1] + c2 * k2[1] + c3 * k3[1] + c4 * k4[1] + c5 * k5[1] + c6 * k6[1]
  I_next = I + c1 * k1[2] + c2 * k2[2] + c3 * k3[2] + c4 * k4[2] + c5 * k5[2] + c6 * k6[2]
  
  if (model_type == "SIR") {
    R_next = N - S_next - I_next
    return(c(S_next, I_next, R_next))
  } else if (model_type == "SIS") {
    return(c(S_next, I_next))
  }
}

# Fonction RK6 avec pas de temps adaptatif pour les modèles SIR et SIS
RK6_adaptive <- function(model_type, beta, gamma, S0, I0, R0, dt_initial, tolerance, max_steps) {
  N <- S0 + I0 + R0
  S <- numeric(max_steps)
  I <- numeric(max_steps)
  R <- NULL
  dt <- dt_initial
  
  S[1] <- S0
  I[1] <- I0
  if (model_type == "SIR") {
    R[1] <- R0
  }
  
  for (step in 1:(max_steps - 1)) {
    repeat {
      current_step <- c(S[step], I[step])
      if (model_type == "SIR") {
        current_step <- c(S[step], I[step], R[step])
      }
      next_step <- RK6_step(model_type, current_step[1], current_step[2], beta, gamma, N, dt)
      half_step1 <- RK6_step(model_type, current_step[1], current_step[2], beta, gamma, N, dt / 2)
      half_step2 <- RK6_step(model_type, half_step1[1], half_step1[2], beta, gamma, N, dt / 2)
      
      error <- estimate_error(next_step, half_step2)
      #nous avons choisie une methode adaptive qui divise le pas sur deux chaque fois que l erreure est sup a la tolerance desire
      if (error > tolerance) {
        dt <- dt / 2
      } else {
        S[step + 1] <- half_step2[1]
        I[step + 1] <- half_step2[2]
        if (model_type == "SIR") {
          R[step + 1] <- half_step2[3]
        }
        if (error < tolerance / 2 && 2 * dt <= dt_initial) {
          dt <- dt * 2
        }
        break
      }
    }
  }
  
  if (model_type == "SIR") {
    return(list(S = S, I = I, R = R))
  } else if (model_type == "SIS") {
    return(list(S = S, I = I))
  }
}

# Paramètres de la simulation pour le modèle SIR
beta_sir <- 0.3
gamma_sir <- 0.1
S0_sir <- 1000
I0_sir <- 1
R0_sir <- 0
dt_initial_sir <- 1
tolerance_sir <- 1
max_steps_sir <- 100

# Exécution de la simulation pour le modèle SIR
result_sir <- RK6_adaptive("SIR", beta_sir, gamma_sir, S0_sir, I0_sir, R0_sir, dt_initial_sir, tolerance_sir, max_steps_sir)

# Paramètres de la simulation pour le modèle SIS
beta_sis <- 0.3
gamma_sis <- 0.1
S0_sis <- 1000
I0_sis <- 1
dt_initial_sis <- 1
tolerance_sis <- 1
max_steps_sis <- 100

# Exécution de la simulation pour le modèle SIS
result_sis <- RK6_adaptive("SIS", beta_sis, gamma_sis, S0_sis, I0_sis, 0, dt_initial_sis, tolerance_sis, max_steps_sis)

# Visualisation des résultats pour le modèle SIR
windows()  # Open a new window
plot(result_sir$S, type='l', col='blue', ylim=c(0, max(c(result_sir$S, result_sir$I, result_sir$R))), ylab='Population', xlab='Time')
lines(result_sir$I, type='l', col='red')
lines(result_sir$R, type='l', col='green')
legend('topright', legend=c('Susceptibles', 'Infectés', 'Rétablis'), col=c('blue', 'red', 'green'), lty=1)
title(main='Simulation du modèle SIR avec méthode RK6 et pas de temps adaptatif')

# Visualisation des résultats pour le modèle SIS
windows()  # Open another new window
plot(result_sis$S, type='l', col='blue', ylim=c(0, max(c(result_sis$S, result_sis$I))), ylab='Susceptible', xlab='Time')
lines(result_sis$I, type='l', col='red')
legend('topright', legend=c('Susceptibles', 'Infectés'), col=c('blue', 'red'), lty=1)
title(main='Simulation du modèle SIS avec méthode RK6 et pas de temps adaptatif')
#NOTE: LES DEUX GRAPHES S'AFFICHE L'UN AU DESSUS D'AUTRE.
