#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   Project: NEXXUS
#   Purpose: Matching Process
#   Version: 0.1
#   Programmer: Xin Huang
#   Date: 09/09/2015
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# read in the fake data for testing

match_raw_data <- read.csv("match_raw_data.csv", header = TRUE, 
                           stringsAsFactors = FALSE)
colnames(match_raw_data) <- tolower(colnames(match_raw_data))




match <- function(data, group, id, mvars, wts, dmaxk, dmax, dist = 1,
                  ncontls = 1, time, transf = 0, seedca, seedco,
                  print = TRUE, rsp_var, rsp_ind, pscore_cut = 0.0001,
                  out, outnmca, outnmcp){
  # @@parameters:
    # DATA    = , /* SAS data set containing cases and potential controls */
    # GROUP   = , /* SAS variable defining cases. Group=1 if case, 0 if control*/
    # ID      = , /* SAS CHARACTER ID variable for the cases and controls*/
    # MVARS   = , /* List of numeric matching variables common to both case and control*/
    # WTS     = , /* List of non-negative weights corresponding to each matching vars*/
    # DMAXK   = , /* List of non-negative values corresponding to each matching vars*/
    # DMAX    = , /* Largest value of Distance(Dij) to be considered as a valid match*/
    # DIST    = 1,/* Type of distance to use 1=absolute difference; 2=Euclidean*/
    # NCONTLS = 1,/* The number of controls to match to each case*/
    # TIME    = , /* Time variable used for risk set matching, only if ctrl time > case being matched*/
    # TRANSF  = 0,/* Whether all matching vars are to be transformed  (0=no, 1=standardize, 2=use ranks)*/
    # SEEDCA  = , /* Seed value used to randomly sort the cases prior to match*/
    # SEEDCO  = , /* Seed value used to randomly sort the controls prior to match*/
    # PRINT   = y,/* Option to print data summary for matched cases */
    # Rsp_var = , /* To include propensity score in match. Response variable to calculate propensity score */
    # Rsp_Ind = , /* Independent variable list to be used for propensity score calculation, dafault=&MVars */
    # Pscore_cut = 0.0001, /*set propensity score diff < 0.005 as valid match */
    # OUT     = __OUT, /* matched data in paired format, &out._Matched in original data format  */
    # OUTNMCA = __NMCA,
    # OUTNMCO = __NMCO);
  
  data = match_raw_data
  group = "case"
  id = "id"
  mvars = c("sex", "age")
  wts = c(1, 1)
  dmaxk = c(0, 2)
  #dmax = 
  dist = 1
  ncontls = 2
  #time
  transf = 0
  seedca = 234
  seedco = 489
  print = TRUE
  rsp_var = "resp"
  rsp_ind = c("sex", "age")
  pscore_cut = 0.0001
  out = "outd"
  outnmca = "matched"
  #outnmcp
  
  # do the logistic regression to get the propensity
  if (!missing(rsp_var)) {
    if (missing(rsp_ind)) {
      rsp_ind <- mvars
    }
    formu <- paste(rsp_var, " ~ ", paste(rsp_ind, collapse = " + "))
    glm_model <- glm(formu, data, family = "binomial")
    p_score <- predict(glm_model, data, type = "response")
    .ndata <- cbind(data, pscore = p_score)
    
    data <- .ndata
    mvars <- c(mvars, "pscore")
    wts <- c(wts, 1)
    if (!missing(dmaxk)) dmaxk <- c(dmaxk, pscore_cut)
  }
  
  bad <- 0
  if (missing(data)) {
    bad <- 1
    stop("ERROR: NO DATASET SUPPLIED") 
  }
  
  if (missing(id)) {
    bad <- 1
    stop("ERROR: NO ID VARIABLE SUPPLIED")  
  }
  
  if (missing(group)) {
    bad <- 1
    stop("ERROR: NO CASE(1)/CONTROL(0) GROUP VARIABLE SUPPLIED")
  }
  
  if (missing(wts)) {
    bad <- 1  
    stop("ERROR: NO WEIGHTS SUPPLIED")
  }
  
  nvar <- length(mvars)
  nwts <- length(wts)
  
  if (nvar != nwts) {
    bad <- 1
    stop("ERROR: #VARS MUST EQUAL #WTS")
  }
  
  nk <- length(dmaxk)
  if (nk > nvar) nk <- nvar
  
  
  v <- mvars
  w <- nwts
  
  if (any(w < 0)) {
    bad = 1
    stop("EERROR: WEIGHTS MUST BE NON-NEGATIVE")
  }
  
  if (nk > 0) {
    k <- dmaxk
    if(any(k < 0)){
      bad = 1
      stop("ERROR: DMAXK VALUES MUST BE NON-NEGATIVE")
    }
  }
  
  ### for match #####
  
  # remove the rows with missing values
  .check <- data
  .check[ ".id"] <- data[, "id"]
  .check <- .check[complete.cases(.check), ]
  
  # standardize the vars
  
  if (transf == 1) {
    .stdzd <- scale(.check[, mvars], center = TRUE, scale = TRUE)
    .caco <- cbind(.check[, setdiff(colnames(.check), mvars)], .stdzd)
  } else if (transf == 0) {
    
  } else {
    .caco <- .check
  }
  
  # for case dataset
  .case <- .caco[.caco[, group] == 1, ]
  
  .case[, ".idca"] <- .case[, id]
  if (!missing(time)) .case[, ".catime"] <- .case[, time]
  
  tmp <- paste(".ca", 1:nvar, sep = "")
  .case[, tmp] <- .case[, v]
  
  if (!missing(seedca)) {
    set.seed(seedca)
    .case[, ".r"] <- rnorm(nrow(.case))
  } else {
    .case[, ".r"] <- 1
  }
  
  if (missing(time)) {
    .case <- .case[, c('.idca', tmp, ".r", mvars)]
  } else {
    .case <- .case[, c('.idca', tmp, ".r", mvars, ".catime")]
  }
  
  nca <- nrow(.case)
  
  # for control dataset
  
  .cont <- .caco[.caco[, group] == 0, ]
  
  .cont[, ".idco"] <- .cont[, id]
  if (!missing(time)) .cont[, ".cotime"] <- .cont[, time]
  
  tmp1 <- paste(".co", 1:nvar, sep = "")
  .cont[, tmp1] <- .cont[, v]
  
  if (!missing(seedco)) {
    set.seed(seedco)
    .cont[, ".r"] <- rnorm(nrow(.cont))
  } else {
    .cont[, ".r"] <- 1
  }
  
  if (missing(time)) {
    .cont <- .cont[, c('.idco', tmp1, ".r", mvars)]
  } else {
    .cont <- .cont[, c('.idco', tmp1, ".r", mvars, ".cotime")]
  }
  
  nco <- nrow(.cont)
  
  bad2 <- 0
  
  if (nco < nca * ncontls) {
    bad2 <- 1
    stop("ERROR: NOT ENOUGH CONTROLS TO MAKE REQUESTED MATCHES")
  }
  
  .cont <- .cont[order(.cont[, c(".r", ".idco")]), ]
  
  #### do the matching #####
  
  # set some flags
  
  .used <- vector("integer", nco)
  .small <- NULL
  .match <- NULL
  
  # distance matrix
  # dm <- as.matrix(.case[, mvars]) %*% as.matrix(t(.cont[, mvars]))
  # tmp <- NULL
  d.matrix <- as.vector(apply(.case[, tmp], 1, 
                              function(x) apply(.cont[, tmp1], 1, 
                                                function(y) x - y)))
  d.matrix <- matrix(d.matrix, ncol = length(tmp), byrow = TRUE)
  d.matrix_m <- apply(d.matrix, 1, function(x) sum(x ^ 2))
  d.matrix_m <- matrix(d.matrix_m, ncol = nrow(.cont), byrow = TRUE)
  
  

  
  
  
  
  
  
 
  
  
  
  
  
  
  
  
  
  
    
  
}

