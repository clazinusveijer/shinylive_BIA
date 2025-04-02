run_simulation_st <- function(c_trt){
  
  #if (!require('pacman')) install.packages('pacman'); library(pacman) 
  #p_load_gh("DARTH-git/darthtools") 
  #p_load("devtools", "dplyr", "scales", "ellipse", "ggplot2", "lazyeval", "igraph", "truncnorm", "ggraph", "reshape2", "stringr", "dampack")
  
  input_parameters <- function(c_trt, modifier_mild_moderate){
    cycle_length   <- 1/365  # cycle length equal to one year
    n_cycles       <- 60    # time horizon, number of cycles
    v_names_states <- c("MV_mild", "ICU_mild", "GW_mild", "REC_mild",   "D_mild","MV_moderate","ICU_moderate", "GW_moderate",  "REC_moderate", "D_moderate","MV_severe","ICU_severe", "GW_severe", "REC_severe", "D_severe")
    n_states       <- length(v_names_states) # number of health states 
    v_names_str   <- c("Standard of care",   # store the strategy names
                       "FX06 treatment") 
    n_str         <- length(v_names_str)     # number of strategies
    
    df_prob_mv_discont <- read.csv("https://raw.githubusercontent.com/clazinusveijer/datasets_CEA_IXION/refs/heads/main/df_prob_mv_discont.csv", sep = ',')
    n_i            <- 2377
    p_i_mild       <- 0.300
    p_i_moderate   <- 0.465
    p_i_severe     <- 0.234
    n_i_mild       <- round(n_i * p_i_mild,0)
    n_i_moderate   <- round(n_i * p_i_moderate,0)
    n_i_severe     <- n_i - n_i_mild - n_i_moderate #round(n_i * p_i_severe,0)
    
    # initial number of patients per health state
    r_i_MV  <- 0.856
    r_i_ICU <- 1-r_i_MV # 0.144 
    n_i_mild_MV  <- round(n_i_mild*r_i_MV,0)
    n_i_mild_ICU <- round(n_i_mild*r_i_ICU,0)
    n_i_moderate_MV  <- round(n_i_moderate*r_i_MV,0)
    n_i_moderate_ICU <- round(n_i_moderate*r_i_ICU,0)
    n_i_severe_MV  <- round(n_i_severe*r_i_MV,0)
    n_i_severe_ICU <- round(n_i_severe*r_i_ICU,0)
    
    modifier_mild_moderate  <- 0.1 #modifier_mild_moderate
    modifier_moderate_severe<- 0.1 #modifier_moderate_severe
    modifier_mild_severe    <- 0.1 #modifier_mild_severe
    
    modifier_MV_ICU_mild      <- 0.1 #modifier_MV_ICU_mild
    modifier_MV_ICU_moderate  <- 0.1 #modifier_MV_ICU_moderate
    modifier_MV_ICU_severe    <- 0.1 #modifier_MV_ICU_severe
    
    # Disease progression up to 7 days post-diagnosis of ARDS
    d_progression_period    <- 7
    r_mild_moderate         <- 0.258
    r_moderate_severe       <- 0.127
    r_mild_severe           <- 0.045
    p_mild_moderate_soc     <- 1-((1-r_mild_moderate)^(1/d_progression_period))
    p_moderate_severe_soc   <- 1-((1-r_moderate_severe)^(1/d_progression_period))
    p_mild_severe_soc       <- 1-((1-r_mild_severe)^(1/d_progression_period))
    
    p_mild_moderate_soc     <- c(p_mild_moderate_soc, 0)
    p_moderate_severe_soc   <- c(p_moderate_severe_soc, 0)
    p_mild_severe_soc       <- c(p_mild_severe_soc, 0)
    
    r_mild_moderate_trt     <- if_else(r_mild_moderate - modifier_mild_moderate < 0, 0, r_mild_moderate - modifier_mild_moderate)
    r_moderate_severe_trt   <- if_else(r_moderate_severe - modifier_moderate_severe < 0, 0, r_moderate_severe - modifier_moderate_severe)
    r_mild_severe_trt       <- if_else(r_mild_severe - modifier_mild_severe < 0, 0, r_mild_severe - modifier_mild_severe)
    
    p_mild_moderate_trt     <- 1-((1-r_mild_moderate_trt)^(1/d_progression_period))
    p_moderate_severe_trt   <- 1-((1-r_moderate_severe_trt)^(1/d_progression_period))
    p_mild_severe_trt       <- 1-((1-r_mild_severe_trt)^(1/d_progression_period))
    
    p_mild_moderate_trt     <- c(p_mild_moderate_trt, 0)
    p_moderate_severe_trt   <- c(p_moderate_severe_trt, 0)
    p_mild_severe_trt       <- c(p_mild_severe_trt, 0)
    
    # 28-day mortality
    r_mort_28_mild          <- 0.296
    r_mort_28_moderate      <- 0.352
    r_mort_28_severe        <- 0.409
    p_D_mild                <- 1 - (1 - r_mort_28_mild)^(1 / 28)
    p_D_moderate            <- 1 - (1 - r_mort_28_moderate)^(1 / 28)
    p_D_severe              <- 1 - (1 - r_mort_28_severe)^(1 / 28)
    
    ### MV LOS
    median_MV_LOS_mild      <- 6
    median_MV_LOS_moderate  <- 8
    median_MV_LOS_severe    <- 11
    
    ### time-dependent MV LOS
    p_MV_ICU_mild_soc         <- df_prob_mv_discont$probability[which(df_prob_mv_discont$severity_group == 'mild')]
    p_MV_ICU_moderate_soc     <- df_prob_mv_discont$probability[which(df_prob_mv_discont$severity_group == 'moderate')]
    p_MV_ICU_severe_soc       <- df_prob_mv_discont$probability[which(df_prob_mv_discont$severity_group == 'severe')]
    
    p_MV_ICU_mild_trt         <- df_prob_mv_discont$probability[which(df_prob_mv_discont$severity_group == 'mild')]*(1+modifier_MV_ICU_mild)
    p_MV_ICU_severe_trt       <- df_prob_mv_discont$probability[which(df_prob_mv_discont$severity_group == 'severe')]*(1+modifier_MV_ICU_moderate)
    p_MV_ICU_moderate_trt     <- df_prob_mv_discont$probability[which(df_prob_mv_discont$severity_group == 'moderate')]*(1+modifier_MV_ICU_severe)
    
    # ICU LOS
    median_ICU_LOS_mild     <- 10
    median_ICU_LOS_moderate <- 12
    median_ICU_LOS_severe   <- 14
    median_ICU_LOS_mild_net     <- median_ICU_LOS_mild - median_MV_LOS_mild
    median_ICU_LOS_moderate_net <- median_ICU_LOS_moderate - median_MV_LOS_moderate
    median_ICU_LOS_severe_net   <- median_ICU_LOS_severe - median_MV_LOS_severe
    
    r_ICU_mild              <- -log(0.5)/median_ICU_LOS_mild_net
    r_ICU_moderate          <- -log(0.5)/median_ICU_LOS_moderate_net
    r_ICU_severe            <- -log(0.5)/median_ICU_LOS_severe_net
    
    p_ICU_GW_mild           <- 1 - exp(-r_ICU_mild)
    p_ICU_GW_moderate       <- 1 - exp(-r_ICU_moderate)
    p_ICU_GW_severe         <- 1 - exp(-r_ICU_severe)
    
    # GW LOS
    median_GW_LOS_mild      <- 23
    median_GW_LOS_moderate  <- 22
    median_GW_LOS_severe    <- 26
    median_GW_LOS_mild_net      <- median_GW_LOS_mild - median_ICU_LOS_mild
    median_GW_LOS_moderate_net  <- median_GW_LOS_moderate - median_ICU_LOS_moderate
    median_GW_LOS_severe_net    <- median_GW_LOS_severe - median_ICU_LOS_severe
    
    r_GW_mild               <- -log(0.5)/median_GW_LOS_mild_net
    r_GW_moderate           <- -log(0.5)/median_GW_LOS_moderate_net
    r_GW_severe             <- -log(0.5)/median_GW_LOS_severe_net
    
    p_GW_REC_mild           <- 1 - exp(-r_GW_mild)
    p_GW_REC_moderate       <- 1 - exp(-r_GW_moderate)
    p_GW_REC_severe         <- 1 - exp(-r_GW_severe)
    
    #### Costs 
    c_MV      <- 2225 # Each additional day on mechanically ventilated ICU patients
    c_ICU     <- 930 # Each additional day on non-mechanically ventilated ICU patients
    c_GW      <- 420  # Each additional day on the general ward for non-ICU patients
    c_REC     <- 0
    c_D       <- 0     # annual cost of being dead
    c_trt     <- c_trt # five-day treatment 600/5
    df_c      <- data.frame(type = c("c_MV", "c_ICU", "c_GW", "c_REC", "c_D", "c_trt"), cost = c(c_MV, c_ICU, c_GW, c_REC, c_D, c_trt))
    
    # sample from age distribution an initial age for every individual
    mean_age <- 61.5
    sd_age <- 14.9
    v_age_init <- rnorm(n_i, mean = mean_age, sd = sd_age)
    min_age <- 16
    max_age <- 99
    df_age <- data.frame(age = round(v_age_init, 0))
    df_age <- df_age %>% mutate(age = if_else(age<min_age, min_age, if_else(age>max_age, max_age, age)))
    df_age <- df_age %>% group_by(age) %>% summarize(prop = n()/n_i)
    v_age_init  <- sample(x = df_age$age, prob = df_age$prop, size = n_i, replace = TRUE) 
    
    v_M_init_mild <- rep(c("MV_mild", "ICU_mild"), times = c(n_i_mild_MV, n_i_mild_ICU))
    df_X_mild <- data.frame(ID = 1:n_i_mild, Severity = "mild", M_init = v_M_init_mild, n_cycles_MV = if_else(v_M_init_mild == "MV_mild", 1, 0),   n_cycles_ICU = 1, trt_YN = 1, mild_moderate_YN = 1, mild_severe_YN = 1, moderate_severe_YN = 2)
    v_M_init_moderate <- rep(c("MV_moderate", "ICU_moderate"), times = c(n_i_moderate_MV, n_i_moderate_ICU))
    df_X_moderate <- data.frame(ID = (n_i_mild+1):(n_i_mild+n_i_moderate), Severity = "moderate", M_init = v_M_init_moderate, n_cycles_MV =   if_else(v_M_init_moderate == "MV_moderate", 1, 0), n_cycles_ICU = 1, trt_YN = 1, mild_moderate_YN = 2, mild_severe_YN = 2, moderate_severe_YN =   1)
    v_M_init_severe <- rep(c("MV_severe", "ICU_severe"), times = c(n_i_severe_MV, n_i_severe_ICU))
    df_X_severe <- data.frame(ID = (n_i_mild+n_i_moderate+1):n_i, Severity = "severe", M_init = v_M_init_severe, n_cycles_MV =   if_else(v_M_init_severe == "MV_severe", 1, 0), n_cycles_ICU = 1, trt_YN = 1, mild_moderate_YN = 2, mild_severe_YN = 2, moderate_severe_YN = 2)
    
    df_X <- rbind(df_X_mild, df_X_moderate, df_X_severe)
    df_X <- df_X %>% mutate(Age = v_age_init,
                            age_group_id = if_else(Age < 65, 1, 2))
    
    inputs <- list(cycle_length = cycle_length
                   , n_cycles = n_cycles
                   , v_names_states = v_names_states
                   , v_names_str = v_names_str
                   , n_str = n_str
                   , n_states = n_states
                   , n_i = n_i
                   , df_c = df_c
                   , df_X = df_X
                   , p_mild_moderate_soc = p_mild_moderate_soc
                   , p_moderate_severe_soc = p_moderate_severe_soc
                   , p_mild_severe_soc = p_mild_severe_soc
                   , p_mild_moderate_trt = p_mild_moderate_trt
                   , p_moderate_severe_trt = p_moderate_severe_trt
                   , p_mild_severe_trt = p_mild_severe_trt
                   , p_D_mild = p_D_mild
                   , p_D_moderate = p_D_moderate
                   , p_D_severe = p_D_severe
                   , p_MV_ICU_mild_soc = p_MV_ICU_mild_soc
                   , p_MV_ICU_moderate_soc = p_MV_ICU_moderate_soc
                   , p_MV_ICU_severe_soc = p_MV_ICU_severe_soc
                   , p_MV_ICU_mild_trt = p_MV_ICU_mild_trt
                   , p_MV_ICU_moderate_trt = p_MV_ICU_moderate_trt
                   , p_MV_ICU_severe_trt = p_MV_ICU_severe_trt
                   , p_ICU_GW_mild = p_ICU_GW_mild
                   , p_ICU_GW_moderate = p_ICU_GW_moderate
                   , p_ICU_GW_severe = p_ICU_GW_severe
                   , p_GW_REC_mild = p_GW_REC_mild
                   , p_GW_REC_moderate = p_GW_REC_moderate
                   , p_GW_REC_severe = p_GW_REC_severe
    )
    return(inputs)
  }
  
  samplev <- function(m_Probs, m = 1) {
    lev <- dimnames(m_Probs)[[2]]  # extract the names of the health states considered for sampling
    n_samp <- nrow(m_Probs)
    u <- runif(n_samp, min = 0, max = 1)
    v_sum_p <- matrixStats::rowCumsums(m_Probs)
    v_cat <- lev[max.col(v_sum_p >= u, ties.method = "first")]
    return(v_cat)
  }
  
  Probs <- function(M_t, t, Trt, input_params, df_X) { 
    #df_X <- as.data.frame(input_params$df_X)
    # Treatment specific transition probabilities
    if (Trt == "Standard of care") {
      p_mild_moderate <- input_params$p_mild_moderate_soc
    } else if (Trt == "FX06 treatment") {
      p_mild_moderate <- input_params$p_mild_moderate_trt
    }  else { 
      warning("Invalid treatment type (mild_moderate)") 
    }
    if (Trt == "Standard of care") {
      p_moderate_severe <- input_params$p_moderate_severe_soc
    } else if (Trt == "FX06 treatment") {
      p_moderate_severe <- input_params$p_moderate_severe_trt
    }  else {
      warning("Invalid treatment type (moderate_severe)") 
    }
    if (Trt == "Standard of care") {
      p_mild_severe <- input_params$p_mild_severe_soc
    } else if (Trt == "FX06 treatment") {
      p_mild_severe <- input_params$p_mild_severe_trt
    }  else {
      warning("Invalid treatment type (mild_severe)") 
    }
    if (Trt == "Standard of care") {
      p_MV_ICU_mild <- input_params$p_MV_ICU_mild_soc
    } else if (Trt == "FX06 treatment") {
      p_MV_ICU_mild <- input_params$p_MV_ICU_mild_trt
    } else {
      warning("Invalid treatment type (MV_ICU_mild)") 
    }
    if (Trt == "Standard of care") {
      p_MV_ICU_moderate <- input_params$p_MV_ICU_moderate_soc
    } else if (Trt == "FX06 treatment") {
      p_MV_ICU_moderate <- input_params$p_MV_ICU_moderate_trt
    } else { 
      warning("Invalid treatment type (MV_ICU_moderate") 
    }
    if (Trt == "Standard of care") {
      p_MV_ICU_severe <- input_params$p_MV_ICU_severe_soc
    } else if (Trt == "FX06 treatment") {
      p_MV_ICU_severe <- input_params$p_MV_ICU_severe_trt
    } else {
      warning("Invalid treatment type (MV_ICU_severe)") 
    }
    
    # create matrix of state transition probabilities
    m_p_t           <- matrix(0, nrow = input_params$n_states, ncol = input_params$n_i)  
    # give the state names to the rows
    rownames(m_p_t) <- input_params$v_names_states                               
    
    ######################## mild
    # transition probabilities when hospitalised in ICU on MV
    m_p_t["MV_mild",  M_t == "MV_mild"]  <- (1 - input_params$p_D_mild - max(p_MV_ICU_mild[df_X$n_cycles_MV]) - max(p_mild_moderate[df_X$mild_moderate_YN]) - max(p_mild_severe[df_X$mild_severe_YN]))
    m_p_t["ICU_mild", M_t == "MV_mild"]  <- max(p_MV_ICU_mild[df_X$n_cycles_MV])
    m_p_t["GW_mild",  M_t == "MV_mild"]  <- 0
    m_p_t["REC_mild", M_t == "MV_mild"]  <- 0
    m_p_t["D_mild",   M_t == "MV_mild"]  <- input_params$p_D_mild    
    
    m_p_t["MV_moderate", M_t == "MV_mild"] <- max(p_mild_moderate[df_X$mild_moderate_YN])
    m_p_t["MV_severe", M_t == "MV_mild"] <- max(p_mild_severe[df_X$mild_severe_YN])
    
    # transition probabilities when hospitalised in ICU without MV
    m_p_t["MV_mild",  M_t == "ICU_mild"] <- 0
    m_p_t["ICU_mild", M_t == "ICU_mild"] <- 1 - input_params$p_D_mild - input_params$p_ICU_GW_mild - max(p_mild_moderate[df_X$mild_moderate_YN]) - max(p_mild_severe[df_X$mild_severe_YN])
    m_p_t["GW_mild",  M_t == "ICU_mild"] <- input_params$p_ICU_GW_mild        
    m_p_t["REC_mild", M_t == "ICU_mild"] <- 0
    m_p_t["D_mild",   M_t == "ICU_mild"] <- input_params$p_D_mild    
    
    m_p_t["ICU_moderate", M_t == "ICU_mild"] <- max(p_mild_moderate[df_X$mild_moderate_YN])
    m_p_t["ICU_severe", M_t == "ICU_mild"] <- max(p_mild_severe[df_X$mild_severe_YN])
    
    # transition probabilities when hospitalised outside ICU
    m_p_t["MV_mild",  M_t == "GW_mild"] <- 0
    m_p_t["ICU_mild", M_t == "GW_mild"] <- 0
    m_p_t["GW_mild",  M_t == "GW_mild"] <- (1 - input_params$p_D_mild - input_params$p_GW_REC_mild)
    m_p_t["REC_mild", M_t == "GW_mild"] <- input_params$p_GW_REC_mild
    m_p_t["D_mild",   M_t == "GW_mild"] <- input_params$p_D_mild    
    
    # transition probabilities when recovered
    m_p_t["MV_mild",  M_t == "REC_mild"] <- 0
    m_p_t["ICU_mild", M_t == "REC_mild"] <- 0
    m_p_t["GW_mild",  M_t == "REC_mild"] <- 0
    m_p_t["REC_mild", M_t == "REC_mild"] <- 1 - input_params$p_D_mild
    m_p_t["D_mild",   M_t == "REC_mild"] <- input_params$p_D_mild  
    
    # transition probabilities when Dead
    m_p_t["MV_mild",  M_t == "D_mild"]  <- 0
    m_p_t["ICU_mild", M_t == "D_mild"]  <- 0
    m_p_t["GW_mild",  M_t == "D_mild"]  <- 0
    m_p_t["REC_mild", M_t == "D_mild"]  <- 0 
    m_p_t["D_mild",   M_t == "D_mild"]  <- 1  
    
    ######################## moderate
    # transition probabilities when hospitalised in ICU on MV
    m_p_t["MV_moderate",  M_t == "MV_moderate"]  <- (1 - input_params$p_D_moderate - max(p_MV_ICU_moderate[df_X$n_cycles_MV]) - max(p_moderate_severe[df_X$moderate_severe_YN]))
    m_p_t["ICU_moderate", M_t == "MV_moderate"]  <- max(p_MV_ICU_moderate[df_X$n_cycles_MV])
    m_p_t["GW_moderate",  M_t == "MV_moderate"]  <- 0
    m_p_t["REC_moderate", M_t == "MV_moderate"]  <- 0
    m_p_t["D_moderate",   M_t == "MV_moderate"]  <- input_params$p_D_moderate 
    
    m_p_t["MV_severe", M_t == "MV_moderate"] <- max(p_moderate_severe[df_X$moderate_severe_YN])
    
    # transition probabilities when hospitalised in ICU without MV
    m_p_t["MV_moderate",  M_t == "ICU_moderate"] <- 0
    m_p_t["ICU_moderate", M_t == "ICU_moderate"] <- (1 - input_params$p_D_moderate - input_params$p_ICU_GW_moderate - max(p_moderate_severe[df_X$moderate_severe_YN]))
    m_p_t["GW_moderate",  M_t == "ICU_moderate"] <- input_params$p_ICU_GW_moderate 
    m_p_t["REC_moderate", M_t == "ICU_moderate"] <- 0
    m_p_t["D_moderate",   M_t == "ICU_moderate"] <- input_params$p_D_moderate    
    
    m_p_t["ICU_severe", M_t == "ICU_moderate"] <- max(p_moderate_severe[df_X$moderate_severe_YN])
    
    # transition probabilities when hospitalised outside ICU
    m_p_t["MV_moderate",  M_t == "GW_moderate"] <- 0
    m_p_t["ICU_moderate", M_t == "GW_moderate"] <- 0
    m_p_t["GW_moderate",  M_t == "GW_moderate"] <- (1 - input_params$p_D_moderate - input_params$p_GW_REC_moderate)
    m_p_t["REC_moderate", M_t == "GW_moderate"] <- input_params$p_GW_REC_moderate
    m_p_t["D_moderate",   M_t == "GW_moderate"] <- input_params$p_D_moderate
    
    # transition probabilities when recovered
    m_p_t["MV_moderate",  M_t == "REC_moderate"] <- 0
    m_p_t["ICU_moderate", M_t == "REC_moderate"] <- 0
    m_p_t["GW_moderate",  M_t == "REC_moderate"] <- 0
    m_p_t["REC_moderate", M_t == "REC_moderate"] <- 1 - input_params$p_D_moderate
    m_p_t["D_moderate",   M_t == "REC_moderate"] <- input_params$p_D_moderate  
    
    # transition probabilities when Dead
    m_p_t["MV_moderate",  M_t == "D_moderate"]  <- 0
    m_p_t["ICU_moderate", M_t == "D_moderate"]  <- 0
    m_p_t["GW_moderate",  M_t == "D_moderate"]  <- 0
    m_p_t["REC_moderate", M_t == "D_moderate"]  <- 0 
    m_p_t["D_moderate",   M_t == "D_moderate"]  <- 1
    
    ######################## severe
    # transition probabilities when hospitalised in ICU on MV
    m_p_t["MV_severe",  M_t == "MV_severe"]  <- (1 - input_params$p_D_severe - max(p_MV_ICU_severe[df_X$n_cycles_MV]))
    m_p_t["ICU_severe", M_t == "MV_severe"]  <- max(p_MV_ICU_severe[df_X$n_cycles_MV])
    m_p_t["GW_severe",  M_t == "MV_severe"]  <- 0
    m_p_t["REC_severe", M_t == "MV_severe"]  <- 0
    m_p_t["D_severe",   M_t == "MV_severe"]  <- input_params$p_D_severe    
    
    # transition probabilities when hospitalised in ICU without MV
    m_p_t["MV_severe",  M_t == "ICU_severe"] <- 0
    m_p_t["ICU_severe", M_t == "ICU_severe"] <- (1 - input_params$p_D_severe - input_params$p_ICU_GW_severe) 
    m_p_t["GW_severe",  M_t == "ICU_severe"] <- input_params$p_ICU_GW_severe 
    m_p_t["REC_severe", M_t == "ICU_severe"] <- 0
    m_p_t["D_severe",   M_t == "ICU_severe"] <- input_params$p_D_severe    
    
    # transition probabilities when hospitalised outside ICU
    m_p_t["MV_severe",  M_t == "GW_severe"] <- 0
    m_p_t["ICU_severe", M_t == "GW_severe"] <- 0
    m_p_t["GW_severe",  M_t == "GW_severe"] <- (1 - input_params$p_D_severe - input_params$p_GW_REC_severe) 
    m_p_t["REC_severe", M_t == "GW_severe"] <- input_params$p_GW_REC_severe
    m_p_t["D_severe",   M_t == "GW_severe"] <- input_params$p_D_severe    
    
    # transition probabilities when recovered
    m_p_t["MV_severe",  M_t == "REC_severe"] <- 0
    m_p_t["ICU_severe", M_t == "REC_severe"] <- 0
    m_p_t["GW_severe",  M_t == "REC_severe"] <- 0
    m_p_t["REC_severe", M_t == "REC_severe"] <- 1 - input_params$p_D_severe
    m_p_t["D_severe",   M_t == "REC_severe"] <- input_params$p_D_severe  
    
    # transition probabilities when Dead
    m_p_t["MV_severe",  M_t == "D_severe"]  <- 0
    m_p_t["ICU_severe", M_t == "D_severe"]  <- 0
    m_p_t["GW_severe",  M_t == "D_severe"]  <- 0
    m_p_t["REC_severe", M_t == "D_severe"]  <- 0 
    m_p_t["D_severe",   M_t == "D_severe"]  <- 1
    
    return(t(m_p_t))
  }
  
  Costs <- function (M_t, Trt, input_params, df_X, df_c) {
    # Arguments:
    # M_t: health state occupied at cycle t (character variable)
    # Returns: 
    # costs accrued in this cycle
    # Trt:  treatment
    
    # Treatment specific transition costs
    if (Trt == "Standard of care") {
      c_trt <- 0
    } else if (Trt == "FX06 treatment") {
      c_trt <- df_c$cost[which(df_c$type == 'c_trt')]
    } 
    
    c_t <- c()
    
    c_t[M_t %in% c("MV_mild", "MV_moderate", "MV_severe")] <- if_else(Trt == "FX06 treatment" & df_X$trt_YN == 1, df_c$cost[which(df_c$type == "c_MV")] + df_c$cost[which(df_c$type == "c_trt")], df_c$cost[which(df_c$type == "c_MV")])
    c_t[M_t %in% c("ICU_mild", "ICU_moderate", "ICU_severe")] <- if_else(Trt == "FX06 treatment" & df_X$trt_YN == 1, df_c$cost[which(df_c$type == "c_ICU")] + df_c$cost[which(df_c$type == "c_trt")], df_c$cost[which(df_c$type == "c_ICU")])
    c_t[M_t %in% c("GW_mild", "GW_moderate", "GW_severe")] <- if_else(Trt == "FX06 treatment" & df_X$trt_YN == 1, df_c$cost[which(df_c$type == "c_GW")] + df_c$cost[which(df_c$type == "c_trt")], df_c$cost[which(df_c$type == "c_GW")])
    c_t[M_t %in% c("REC_mild", "REC_moderate", "REC_severe")]  <- df_c$cost[which(df_c$type == "c_REC")] 
    c_t[M_t %in% c("D_mild", "D_moderate", "D_severe")]    <- df_c$cost[which(df_c$type == "c_D")]  
    
    return(c_t)  # return costs accrued this cycle
  }
  
  MicroSim <- function(Trt, seed, input_params, df_X, df_c) {
    
    set.seed(seed)
    n_states <- length(input_params$v_names_states) # the number of health states
    
    m_M <- m_C <-  matrix(nrow = input_params$n_i, ncol = input_params$n_cycles + 1, 
                          dimnames = list(paste("ind"  , 1:input_params$n_i, sep = " "), 
                                          paste("cycle", 0:input_params$n_cycles, sep = " ")))  
    
    m_M [, 1] <- as.character(input_params$df_X$M_init) # initial health state at cycle 0 for individual i
    # calculate costs per individual during cycle 0
    m_C[, 1]  <- Costs(m_M[, 1], Trt=Trt, input_params = input_params, df_X = df_X, df_c = df_c)     
    
    # open a loop for time running cycles 1 to n_cycles 
    for (t in 1:input_params$n_cycles) {
      # calculate the transition probabilities for the cycle based on health state t
      m_P <- Probs(m_M[, t], t, Trt, input_params, df_X)
      # check if transition probabilities are between 0 and 1
      #check_transition_probability(m_P, verbose = F)
      # check if each of the rows of the transition probabilities matrix sum to one
      #check_sum_of_transition_array(m_P, n_rows = input_params$n_i, n_cycles = input_params$n_cycles, verbose = F)
      
      # sample the next health state and store that state in matrix m_M
      m_M[, t + 1]  <- samplev(m_P, 1)    
      # calculate costs per individual during cycle t + 1
      m_C[, t + 1]  <- Costs(m_M[, t + 1], Trt=Trt, input_params = input_params, df_X = df_X, df_c = df_c)     
      
      # update time since illness onset for t + 1 
      df_X$mild_moderate_YN <- df_X$mild_severe_YN <- if_else(t+1 < 7 & df_X$Severity == "mild", 1, 2)
      df_X$moderate_severe_YN <- if_else(t+1 < 7 & df_X$Severity == "moderate", 1, 2)
      
      df_X$n_cycles_MV <- if_else(grepl("MV", m_M[, t + 1]), df_X$n_cycles_MV +1, df_X$n_cycles_MV)
      
      df_X$n_cycles_ICU <- if_else(m_M[, t + 1] %in% c("MV_mild", "MV_moderate", "MV_severe", "ICU_mild", "ICU_moderate", "ICU_severe"), 
                                   df_X$n_cycles_ICU +1,
                                   df_X$n_cycles_ICU)
      
      df_X$trt_YN <- if_else(m_M[, t + 1] %in% c("MV_mild", "MV_moderate", "MV_severe",
                                                 "ICU_mild", "ICU_moderate", "ICU_severe",
                                                 "GW_mild", "GW_moderate", "GW_severe")
                             & t+1 < 5
                             , 1
                             , 0)
      
      # Display simulation progress
      if(t/(input_params$n_cycles/10) == round(t/(input_params$n_cycles/10), 0)) { # display progress every 10%
        cat('\r', paste(t/input_params$n_cycles * 100, "% done", sep = " "))
      }
      
    } # close the loop for the time points 
    
    tc      <- rowSums(m_C)   # total cost per individual
    tc_hat  <- mean(tc)       # average cost 
    # store the results from the simulation in a list
    results <- list(df_X = df_X, m_M = m_M, m_C = m_C, tc = tc, tc_hat = tc_hat)   
    
    return(results)  # return the results
  }
  
  run_simulation <- function(c_trt){
    input_params <- input_parameters(c_trt=c_trt)
    df_X <- as.data.frame(input_params$df_X)
    df_c <- as.data.frame(input_params$df_c)
    outcomes_SoC   <- MicroSim(Trt="Standard of care", seed = 77, input_params = input_params, df_X = df_X, df_c = df_c)
    outcomes_trt   <- MicroSim(Trt="FX06 treatment", seed = 77, input_params = input_params, df_X = df_X, df_c = df_c)
    outcomes <- list(outcomes_SoC = outcomes_SoC, outcomes_trt = outcomes_trt)
    return(outcomes)
  }
  outcomes <- run_simulation(c_trt)
  return(outcomes)
}