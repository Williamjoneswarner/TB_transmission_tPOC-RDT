#Creating the parametres and initial values. 

#' - inf.prob:         Probability of transmission per contact
#' - act.rate:         Number of contacts per person per time step
#' - piu.rate:         Rate of progression from latent to early infectious (E → P)
#' - delay.rate:       Rate of progression from preclinical to clinical TB (P → I)
#' - p_diag_coverage:  Proportion of individuals in the preclinical stage (P) eligible for diagnosis
#' - p_diag_rate:      Probability of diagnosis among eligible individuals in the preclinical stage (P → D)
#' - i_diag_coverage:  Proportion of individuals in the clinical stage (I) eligible for diagnosis
#' - i_diag_rate:      Probability of diagnosis among eligible individuals in the clinical stage (I → D)
#' - treat.rate:       Probability of treatment uptake once diagnosed (D → T)
#' - rec.rate:         Rate of recovery following treatment (T → R)
#' - recid.rate:       Rate of relapse or TB-related death, balanced by replacement into susceptible pool (R → S)
#' 

params <- list(
  inf.prob     = 0.0024,
  act.rate     = 5,#this is taking the average between rural and urban settings.
  piu.rate     = 1 / 365, #1/365 possible
  delay.rate   = 1 / 45,
  p_diag_coverage = 0.5, #proportion of population where the diagnostic acts againstP_ og was 0.1.
  p_diag_rate  = 1 / 120, #og was 1/183 putting up to 1 / 120 to match current diagnostic rate
  i_diag_coverage = 0.9, # og was 0.76 putting up to 0.9
  i_diag_rate  = 1.7 / 120, #no change in diagnosis rate, as it is unlikley to impact this.
  treat.rate   = 1 / 10,
  rec.rate     = 1 / 180,
  recid.rate   = 0.6 / 365
)

# params <- list(
#   inf.prob     = 0.0024,
#   act.rate     = 5,
#   piu.rate     = 1 / 365, #1/365 possible
#   delay.rate   = 1 / 45,
#   p_diag_rate  = 0.2 / 183,
#   i_diag_rate  = 1.2 / 120,
#   treat.rate   = 1 / 10,
#   rec.rate     = 1 / 180,
#   recid.rate   = 0.2 / 365
# )

# Create time column
time_points <- 1:630

# Create initial row values
initial_values <- data.frame(
  s.num = 9000, #suseptable number
  e.num = 1000, #exposed number
  p.num = 1000, #preclinical number
  i.num = 1000, #infectious number
  d.num = 1000, #dioagnosed number
  t.num = 1000, #treated number
  r.num = 1000 #recovered number
)

generate_initial_state_column <- function(initial_values) {
  values <- unlist(initial_values)
  total_n <- sum(values)
  state_labels <- gsub("\\.num", "", names(values))
  
  # Create vector of states according to counts (in order)
  all_states <- unlist(mapply(rep, state_labels, values))
  
  # Create data frame with one column
  df <- data.frame(state = all_states, stringsAsFactors = FALSE)
  
  return(df)
}

state_df <- generate_initial_state_column(initial_values)
head(state_df)
table(state_df)

#1) ##infection module - s > e

infection_module <- function(state_df, params, timepoint) {
  current_col <- if (timepoint == 1) "state" else paste0("timepoint_", timepoint)
  states <- state_df[[current_col]]
  
  # Number susceptible
  n_sus <- sum(states == "s")
  
  # Initialize next column as NA
  new_col <- rep(NA, length(states))
  
  if (n_sus > 0) {
    # Calculate force of infection
    infectious_states <- c("i", "d")
    prop_infectious <- sum(states %in% infectious_states) / length(states)
    lambda <- params$inf.prob * params$act.rate * prop_infectious
    
    # Identify susceptible indices
    sus_indices <- which(states == "s")
    
    # Run Bernoulli trials
    rand_vals <- runif(n_sus)
    infected <- sus_indices[rand_vals < lambda]
    
    # Update states: infected → 'e', others stay 's'
    new_col[sus_indices] <- "s"
    new_col[infected] <- "e"
  }
  
  # Add new timepoint column
  new_col_name <- paste0("timepoint_", timepoint + 1)
  state_df[[new_col_name]] <- new_col
  
  return(state_df)
}

state_df <- infection_module(state_df, params, timepoint = 1)
colnames(state_df)
table(state_df$state)
table(state_df$timepoint_2)

#2) function for transitioning exposed to preclinical. - e > p

progress_E_to_P <- function(state_df, params, timepoint) {
  current_col <- if (timepoint == 1) "state" else paste0("timepoint_", timepoint)
  next_col <- paste0("timepoint_", timepoint + 1)
  
  current_states <- state_df[[current_col]]
  
  # Create next column if not already present
  if (!next_col %in% names(state_df)) {
    state_df[[next_col]] <- rep(NA, length(current_states))
  }
  
  # Identify rows that are 'e' and not yet updated in the next column
  to_update <- which(is.na(state_df[[next_col]]) & current_states == "e")
  n_e <- length(to_update)
  
  if (n_e == 0) return(state_df)
  
  # Progression probability from E to P
  p_progress <- params$piu.rate
  
  # Bernoulli trials
  rand_vals <- runif(n_e)
  progressing <- to_update[rand_vals < p_progress]
  staying <- setdiff(to_update, progressing)
  
  # Update states
  state_df[[next_col]][progressing] <- "p"
  state_df[[next_col]][staying] <- "e"
  
  return(state_df)
}

state_df <- progress_E_to_P(state_df, params, timepoint = 1)
table(state_df$state)
table(state_df$timepoint_2)

#3) Preclinical to Infectious  

progress_P <- function(state_df, params, timepoint) {
  # Determine relevant column names
  current_col <- if (timepoint == 1) "state" else paste0("timepoint_", timepoint)
  next_col <- paste0("timepoint_", timepoint + 1)
  
  current_states <- state_df[[current_col]]
  
  # Create the next column if not present
  if (!next_col %in% names(state_df)) {
    state_df[[next_col]] <- rep(NA, length(current_states))
  }
  
  # Identify "p" individuals who haven't yet updated
  to_update <- which(is.na(state_df[[next_col]]) & current_states == "p")
  n_p <- length(to_update)
  if (n_p == 0) return(state_df)
  
  # Step 1: randomly select a subset of to_update to be eligible for diagnosis
  diag_eligible_n <- round(params$p_diag_coverage * n_p)
  diag_eligible <- sample(to_update, size = diag_eligible_n, replace = FALSE)
  
  # Apply P → D transition to eligible subset
  p_diag <- params$p_diag_rate
  rand_vals_diag <- runif(length(diag_eligible))
  diag_ids <- diag_eligible[rand_vals_diag < p_diag]
  state_df[[next_col]][diag_ids] <- "d"
  
  # Step 2: everyone else (not diagnosed) is evaluated for P → I
  remaining <- setdiff(to_update, diag_ids)
  if (length(remaining) > 0 && params$delay.rate > 0) {
    rand_vals_inf <- runif(length(remaining))
    inf_ids <- remaining[rand_vals_inf < params$delay.rate]
    state_df[[next_col]][inf_ids] <- "i"
    
    # Step 3: Leave remaining as "p"
    still_p <- setdiff(remaining, inf_ids)
    state_df[[next_col]][still_p] <- "p"
  } else {
    state_df[[next_col]][remaining] <- "p"
  }
  
  return(state_df)
}

state_df <- progress_P(state_df, params, timepoint = 1)
table(state_df$state)
table(state_df$timepoint_2)

#4) Transitioning infectious to diagnosed using conventional diagnostics. 

progress_I_to_D <- function(state_df, params, timepoint) {
  # Determine column names for current and next timepoint
  current_col <- if (timepoint == 1) "state" else paste0("timepoint_", timepoint)
  next_col <- paste0("timepoint_", timepoint + 1)
  
  current_states <- state_df[[current_col]]
  
  # Create the next column if not already present
  if (!next_col %in% names(state_df)) {
    state_df[[next_col]] <- rep(NA, length(current_states))
  }
  
  # Identify rows that are "i" and not yet updated
  to_update <- which(is.na(state_df[[next_col]]) & current_states == "i")
  n_i <- length(to_update)
  if (n_i == 0) return(state_df)
  
  # Step 1: randomly select a subset of to_update eligible for diagnosis
  diag_eligible_n <- round(params$i_diag_coverage * n_i)
  diag_eligible <- sample(to_update, size = diag_eligible_n, replace = FALSE)
  
  # Step 2: Apply diagnosis rate among eligible
  i_diag <- params$i_diag_rate
  rand_vals <- runif(length(diag_eligible))
  diagnosed <- diag_eligible[rand_vals < i_diag]
  
  # Step 3: Everyone else stays as "i"
  staying <- setdiff(to_update, diagnosed)
  
  # Update states
  state_df[[next_col]][diagnosed] <- "d"
  state_df[[next_col]][staying] <- "i"
  
  return(state_df)
}

state_df <- progress_I_to_D(state_df, params, timepoint = 1)
table(state_df$state)
table(state_df$timepoint_2)

#4) initiation of treatment on the diagnosed individuals.

progress_D_to_T <- function(state_df, params, timepoint) {
  # Determine current and next timepoint column names
  current_col <- if (timepoint == 1) "state" else paste0("timepoint_", timepoint)
  next_col <- paste0("timepoint_", timepoint + 1)
  
  current_states <- state_df[[current_col]]
  
  # Create next column if it does not exist
  if (!next_col %in% names(state_df)) {
    state_df[[next_col]] <- rep(NA, length(current_states))
  }
  
  # Identify rows currently "d" and not yet updated for next timepoint
  to_update <- which(is.na(state_df[[next_col]]) & current_states == "d")
  n_d <- length(to_update)
  if (n_d == 0) return(state_df)  # no one to update
  
  # Treatment initiation probability
  p_treat <- params$treat.rate
  
  # Bernoulli trials for treatment start
  rand_vals <- runif(n_d)
  treated <- to_update[rand_vals < p_treat]
  staying <- setdiff(to_update, treated)
  
  # Update states
  state_df[[next_col]][treated] <- "t"
  state_df[[next_col]][staying] <- "d"
  
  return(state_df)
}

state_df <- progress_D_to_T(state_df, params, timepoint = 1)
table(state_df$state)
table(state_df$timepoint_2)

#5) treatment success rate.

progress_T_to_R <- function(state_df, params, timepoint) {
  # Determine current and next timepoint columns
  current_col <- if (timepoint == 1) "state" else paste0("timepoint_", timepoint)
  next_col <- paste0("timepoint_", timepoint + 1)
  
  current_states <- state_df[[current_col]]
  
  # Create next column if not present
  if (!next_col %in% names(state_df)) {
    state_df[[next_col]] <- rep(NA, length(current_states))
  }
  
  # Identify rows currently "t" and not updated yet
  to_update <- which(is.na(state_df[[next_col]]) & current_states == "t")
  n_t <- length(to_update)
  if (n_t == 0) return(state_df)
  
  # Recovery probability
  p_recover <- params$rec.rate
  
  # Bernoulli trials for recovery
  rand_vals <- runif(n_t)
  recovered <- to_update[rand_vals < p_recover]
  staying <- setdiff(to_update, recovered)
  
  # Update states
  state_df[[next_col]][recovered] <- "r"
  state_df[[next_col]][staying] <- "t"
  
  return(state_df)
}

state_df <- progress_T_to_R(state_df, params, timepoint = 1)
table(state_df$state)
table(state_df$timepoint_2)

#6) Rec moving to suseptable to represent death rate (treamtent failure) 
# but also relapse rate. 

rec_to_susceptible <- function(state_df, params, timepoint) {
  # Determine current and next timepoint columns
  current_col <- if (timepoint == 1) "state" else paste0("timepoint_", timepoint)
  next_col <- paste0("timepoint_", timepoint + 1)
  
  current_states <- state_df[[current_col]]
  
  # Create next column if not present
  if (!next_col %in% names(state_df)) {
    state_df[[next_col]] <- rep(NA, length(current_states))
  }
  
  # Identify rows currently "r" and not yet updated
  to_update <- which(is.na(state_df[[next_col]]) & current_states == "r")
  n_r <- length(to_update)
  if (n_r == 0) return(state_df)
  
  # Relapse probability
  p_recur <- params$recid.rate
  
  # Bernoulli trials for relapse
  rand_vals <- runif(n_r)
  relapsed <- to_update[rand_vals < p_recur]
  staying <- setdiff(to_update, relapsed)
  
  # Update states
  state_df[[next_col]][relapsed] <- "s"
  state_df[[next_col]][staying] <- "r"
  
  return(state_df)
}

state_df <- rec_to_susceptible(state_df, params, timepoint = 1)
table(state_df$state)
table(state_df$timepoint_2)

### Now lopping functions depending on timepoints.
# infection_module
# progress_E_to_P
# progress_P
# progress_I_to_D
# progress_D_to_T
# progress_T_to_R
# rec_to_susceptible

run_simulation <- function(state_df, params, n_steps) {
  for (t in 1:n_steps) {
    # 1. Infection module (assuming this updates state_df)
    state_df <- infection_module(state_df, params, timepoint = t)
    
    # 2. E to P progression
    state_df <- progress_E_to_P(state_df, params, timepoint = t)
    
    # 3. P progression (P -> I or D)
    state_df <- progress_P(state_df, params, timepoint = t)
    
    # 4. I to D diagnosis
    state_df <- progress_I_to_D(state_df, params, timepoint = t)
    
    # 5. D to T treatment initiation
    state_df <- progress_D_to_T(state_df, params, timepoint = t)
    
    # 6. T to R recovery
    state_df <- progress_T_to_R(state_df, params, timepoint = t)
    
    # 7. R to S relapse
    state_df <- rec_to_susceptible(state_df, params, timepoint = t)
  }
  
  return(state_df)
}

###initial pop values

# Create initial row values

initial_state_df <- data.frame(
  state = c(
    rep("s", 8382),#suseptable number
    rep("e", 457),#exposed number
    rep("p", 54),#preclinical number
    rep("i", 111),#infectious number
    rep("d", 10),#dioagnosed number
    rep("t", 236),#treated number
    rep("r", 750)#recovered number
  )
)

time_points <- 3650*2

# Applying the function. 
final_state <- run_simulation(state_df = initial_state_df, params = params, n_steps = time_points)

# Define the state labels you're interested in
states <- c("s", "e", "p", "i", "d", "t", "r")

# Get the timepoint columns (assumes names like "timepoint_1", ..., "timepoint_630")
names(final_state)[names(final_state) == "state"] <- "timepoint_1"
timepoint_cols <- grep("^timepoint_", names(final_state), value = TRUE)

# Initialize a summary data frame
summary_df <- data.frame(timepoint = seq_along(timepoint_cols))

# For each state, count occurrences per column
for (state in states) {
  summary_df[[state]] <- sapply(timepoint_cols, function(col) sum(final_state[[col]] == state, na.rm = TRUE))
}

library(ggplot2)

ggplot(summary_df, aes(x = timepoint)) +
  geom_line(aes(y = s, color = "Susceptible")) +
  geom_line(aes(y = e, color = "Exposed (Latent)")) +
  geom_line(aes(y = p, color = "Preclinical")) +
  geom_line(aes(y = i, color = "Infectious (ATBI)")) +
  geom_line(aes(y = d, color = "Diagnosed")) +
  geom_line(aes(y = t, color = "Treated")) +
  geom_line(aes(y = r, color = "Recovered")) +
  scale_color_manual(
    name = "Compartment",
    values = c(
      "Susceptible" = "blue",
      "Exposed (Latent)" = "purple",
      "Preclinical" = "orange",
      "Infectious (ATBI)" = "red",
      "Diagnosed" = "brown",
      "Treated" = "pink",
      "Recovered" = "green"
    )
  ) +
  labs(
    title = "Dynamics of compartments over time - Baseline",
    x = "Timepoint (Days)",
    y = "Count"
  )

summary_df[3300*2,]

