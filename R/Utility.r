#**************************************************************************
#
# Copyright Aurobindh Kalathil Puthanpura, ETA, PSU, 2015
# Use granted under BSD license terms
#
# R DEA M u l t ip l ie r Model Package
#
# Author: Aurobindh Kalathil Puthanpura
# Revision: 1
# Date: 2015-20-9
#
# U t i l i t y f i l e
# internal use only
#
#***************************************************************************

#
# Package Global variables
#
#multiplierdea.version = "0.1.0"
#multiplierdea.id = paste0("multiplierdea version ", multiplierdea.version)

# Load ipsolverAPI library

# .loadLibrary <- function() {
#   if (requireNamespace("lpSolveAPI",quietly = TRUE))
#   {
#     library("lpSolveAPI")
#   } else
#   {
#     install.packages("lpSolveAPI",repos = "http://cran.us.r-project.org")
#     library("lpSolveAPI")
#   }
#
#   if (requireNamespace("lpSolveAPI" ,quietly = TRUE))
#   {
#     return(TRUE)
#   } else
#   {
#     return(FALSE)
#   }
#
# }

#
# Options Strings
#

options.orientation.l <- c ("input", "output")
options.rts.l <- c ("vrs" , "crs")
options.phase.l <- c("mal", "ben")

#
# Solver Status Codes
#
dict.solveStatus<-cbind(key= c(0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13),
            desc = c("optimal solution found","the model is sub-optimal", "the model is infeasible",
                     "the model is unbounded", "the model is degenerate", "numerical failure encountered",
                     "process aborted", "timeout", "the model was solved by presolve","the branch and bound routine failed",
                     "the branch and bound was stopped because of a break-at-first or break-at-value",
                     "a feasible branch and bound solution was found","no feasible branch and bound solution was found"))

# Function: .checkDate
# Check for valid types
# Check add names to unnamed variables

.checkData <- function(x, name) {
  # Is value a dataframe, array or vector?
  if (!is.data.frame(x) && !is.array (x) && (length(x) < 2))
    stop(name,
         " is not a matrix, array (1 or 2 dimensions) or data.frame",
         call. =
           FALSE)
  if (length(dim(x)) > 2)
    stop(name, " is greater then 2 dimensions", call. = FALSE)

  # If data.frame - convert to array
  if (is.data.frame(x))
    x <- data.matrix(x)

  # If vector convert to a 1 x N array
  if (is.vector (x))
    X <- array(x, c (length(x), 1))

  # Check that a l l numeric
  if (!is.numeric(x))
    stop(name, " must be numeric", call. = FALSE)

  # Add DMU names i f empty
  if (length(dim(x)) >= 1 && is.null(rownames(x)))
    rownames(x) <- paste0("DMU", 1:nrow(x))

  # Add Col names i f empty
  if (length(dim(x)) >= 2 && is.null(colnames(x)))

    colnames(x) <- paste0(name, 1:ncol(x))
  return(x)
}


# Function .CheckDataGood
# Check input & output data and warn about problems for dea
.checkDataGood <- function(x, y) {
  status <- TRUE

  # Check for any zero's
  # check for na' s *
  # Check for no positive values on input & output
  if (any(x == 0, na.rm = TRUE)) {
    cat(
      "Warning, data has DMU's with inputs that are zero, th is may cause numerical
      problems\n"
    )
    status <- FALSE
  }

  if (any(y == 0, na.rm = TRUE)) {
    cat(
      "Warning, data has DMU's with outputs that are zero, this may cause
      numerical problems\n"
    )
    status <- FALSE
  }

  for (k in (1:nrow (x))) {
    if (all(x[k, ] <= 0, na.rm = TRUE) & all(y[k, ] <= 0, na.rm = TRUE)) {
      cat(
        "warning, dmu # " ,
        k ,
        "has no positive inputs or outputs, this may cause
        numerical problems\n"
      )
      status <- FALSE
    }
  }
  return(status)
  }

# Function .checkWeightRestriction
# Check if the weight restriction dataframe has values
.checkWeightRestriction <- function(x, y, weightRestriction) {
  #Get the coloum names from the IP and OP data set
  colValues <- c(colnames(x), colnames(y))

  #Check lower bound values for NA
  if (any(is.na(weightRestriction[, 1]), na.rm = FALSE))
  {
    stop("Lower bound for the weight restriction table has NA values",
         call. = FALSE)
  }

  #Check upper bound values for NA
  # if(any(is.na(weightRestriction[,4]),na.rm = FALSE))
  # {
  #   stop("Upper bound for the weight restriction table has NA values", call. = FALSE)
  # }


  #Check numerator column names
  for (i in 1:nrow(weightRestriction)) {
    if (length(which(colValues == weightRestriction[i, 2])) == 0)
    {
      stop(
        "Numerator values does not match with the column names from input and output values",
        call. = FALSE
      )
    }
  }

  #Check denominator column names
  for (i in 1:nrow(weightRestriction)) {
    if (length(which(colValues == weightRestriction[i, 3])) == 0)
    {
      stop(
        "Denominator values does not match with the column names from input and output values",
        call. = FALSE
      )
    }
  }

}

# Function .checkoption()
# Check that value in legal l i s t options
# If on legal list, return lowercase short form option
#
.checkoption <- function(value, name, options.l = c(true, false)) {
  if (length(value) != 1) {
    stop("illegal vector for ",
         name,
         " option must be single value",
         call. = false)
  }
  # If options are character
  if (is.character(options.l[1])) {
    if (is.character(value)) {
      tmp.value <- tolower(value)
      # if (tmp.value %in% options.1)
      i <- charmatch(tmp.value, options.l, nomatch = -1)
      if (i > 0)
        return(options.l[i])
    }
  } else if (is.logical(Options.l[1])) {
    # Logical options
    if (is.logical(value)) {
      return(value)
    }
  } else {
    # Numeric options
    if (is.numeric(value) && is.finite(value) && (value >= 0)) {
      return(value)
    }
    options.l <- c("numeric >= 0")
  }
  stop(
    "Illegal value=",
    value,
    " for ",
    name,
    "; legal values are: ",
    paste(options.l, collapse = ", ") ,
    call. = FALSE
  )
}

.Mult_lP <- function (x_t, y_t, rts = 'crs', weightRestriction) {
  Result.rts <- 'crs'
  if (rts == 'vrs') {
    Result.rts <- 'vrs'
  } else
  {
    Result.rts <- 'crs'
  }

  Result.orientation <- "Input"
  Result.Inputvalues <- x_t
  Result.Outputvalues <- y_t
  Result.eff <- matrix(nrow = nrow(x_t))
  Result.Lambda <- matrix(nrow = nrow(x_t), ncol = nrow(x_t))
  Result.vx <- matrix(nrow = nrow(x_t), ncol = ncol(x_t))
  Result.uy <- matrix(nrow = nrow(y_t), ncol = ncol(y_t))
  Result.Free_weights <- matrix(nrow = nrow(x_t), ncol = 2)
  Result.StatusData <- vector(mode="logical", nrow(x_t))



  #free sign variable
  fs <- 0

  switch(rts,
         'crs' = fs <- 0,
         'vrs' = fs <- 2,
         fs <- 2)

  dmu <- row.names(x_t)

  # Assign names to Matrix for efficiency
  rownames(Result.eff) <- dmu
  colnames(Result.eff)[1] <- c("Eff")

  # Assign names to Matrix for efficiency
  rownames(Result.Lambda) <- dmu
  colnames(Result.Lambda) <- dmu

  # Assign names to Matrix for IP Weights
  rownames(Result.vx) <- dmu
  colnames(Result.vx) <- colnames(x_t)

  # Assign names to Matrix for OP Weights
  rownames(Result.uy) <- dmu
  colnames(Result.uy) <- colnames(y_t)

  # Assign names to Matrix for Free Weights
  rownames(Result.Free_weights) <- dmu
  colnames(Result.Free_weights) <- c("F1", "F2")



  for (d in 1:length(dmu)) {
    ip_v <- ncol(x_t)
    op_v <- ncol(y_t)

    # Set objective function for DMU under consideration
    c <- 0
    for (i in 1:ncol(y_t)) {
      c[i] <- c(y_t[d, i])
    }


    # if vrs add the free variables else omit
    if (rts == 'vrs') {
      # Adding the 1 and -1 for free variable
      if (length(c) != 0) {
        c[length(c) + 1] <- 1
        c[length(c) + 1] <- -1
      }
    }

    # indices for objective function
    ind_obj <- 0
    for (i in 1:ncol(y_t)) {
      ind_obj[i] <- c(ip_v + i)
    }

    # if vrs add the free variables indices else omit
    if (rts == 'vrs') {
      # indices for objective function for free variable
      if (length(ind_obj) != 0) {
        ind_obj[length(ind_obj) + 1] <- ind_obj[length(ind_obj)] + 1
        ind_obj[length(ind_obj) + 1] <- ind_obj[length(ind_obj)] + 1
      }
    }

    lps.model <- make.lp(0, ip_v + op_v + fs)
    set.objfn(lps.model, c, indices = ind_obj)

    # Input constraint for the DMU under consideration
    ip_con_dmu <- 0
    for (i in 1:ncol(x_t)) {
      ip_con_dmu[i] <- x_t[d, i]
    }

    #indices for DMU input constraints
    ind_ip_con <- 0
    for (i in 1:ncol(x_t)) {
      ind_ip_con[i] <- i
    }

    add.constraint(lps.model, ip_con_dmu, "=", rhs = 1, ind_ip_con)

    # Multiply input values by -1 for building the constraints for all OMU's
    x_temp <- -1 * x_t

    # Construct constraints for all dmus
    for (i in 1:nrow(x_t)) {
      dmu_con <- 0
      for (j in 1:ncol(x_temp)) {
        dmu_con[j] <- x_temp[i, j]
      }
      for (k in 1:ncol(y_t)) {
        dmu_con[ncol(x_temp) + k] <- y_t[i, k]
      }

      # If rts is vrs then add the free variable else omit
      if (rts == 'vrs') {
        #Adding free variable in the constraint
        if (length(dmu_con) != 0) {
          dmu_con[length(dmu_con) + 1] <- 1
          dmu_con[length(dmu_con) + 1] <- -1
        }

      }
      add.constraint(lps.model, dmu_con, "<=", rhs = 0)
    }
    #Weight Restrictions
    if (!missing(weightRestriction)) {
      colvalues <- c(colnames(x_t), colnames(y_t))

      for (i in 1:nrow(weightRestriction))
      {
        #First constraint
        add.constraint(lps.model,
                       c(weightRestriction[i, 1], -1),
                       ">=",
                       0,
                       indices = c(
                         which(colvalues == weightRestriction[i, 2]),
                         which(colvalues == weightRestriction[i, 3])
                       ))
        #Second constraint
        # Check if upper bound exists (Length ==4)
        if (length(weightRestriction) == 4) {
          #Check if upper bound is not Inf, NA, or NaN

          if (!is.nan(weightRestriction[i, 4]) &&
              !is.na(weightRestriction[i, 4]) &&
              !is.infinite(weightRestriction[i, 4])) {
            add.constraint(
              lps.model,
              c(1, -1 * weightRestriction[i, 4]),
              "<=",
              0,
              indices = c(
                which(colvalues == weightRestriction[i, 2]),
                which(colvalues == weightRestriction[i, 3])
              )
            )
          }
        }
      }
    }

    lp.control(lps.model, sense = 'max')

    #print(lps.model)
    Result.StatusData[d]<-solve(lps.model)

    Result.eff[d] <- get.objective(lps.model)
    temp_lamdba <- get.dual.solution(lps.model)
    Result.Lambda[d, ] <- temp_lamdba[3:(length(dmu) + 2)]


    temp_weight <- get.variables(lps.model)

    #Store the IP & OP weights in respective Matrices
    Result.vx[d, 1:ncol(x_t)] <- temp_weight[1:ncol(x_t)]
    Result.uy[d, 1:ncol(y_t)] <-
      temp_weight[(ncol(x_t) + 1):(ncol(x_t) + ncol(y_t))]
    if (rts == 'vrs') {
      Result.Free_weights[d, 1:2] <-
        temp_weight[(ncol(x_t) + ncol(y_t) + 1):length(temp_weight)]
    }
  }

  # Calculating the transpose of Lambda values for HCU calculation
  lambda_Transpose <- t(Result.Lambda)

  # Calculating HCU data for Input
  Hcu_input <- matrix(nrow = ncol(lambda_Transpose), ncol = ncol(x_t))
  colnames(Hcu_input) <- colnames(x_t)
  rownames(Hcu_input) <- rownames(lambda_Transpose)

  for (i in 1:ncol(lambda_Transpose)) {
    e <- x_t * lambda_Transpose[, i]
    for (j in 1:ncol(e)) {
      Hcu_input[i, j] <- sum(e[, j])
    }
  }

  # Calculating HCU data for Output
  Hcu_output <- matrix(nrow = ncol(lambda_Transpose), ncol = ncol(y_t))
  colnames(Hcu_output) <- colnames(y_t)
  rownames(Hcu_output) <- rownames(lambda_Transpose)

  for (i in 1:ncol(lambda_Transpose)) {
    e <- y_t * lambda_Transpose[, i]
    for (j in 1:ncol(e)) {
      Hcu_output[i, j] <- sum(e[, j])
    }
  }


  #Get the Status code and description from the LP model
  Result.StatusDesc<-sapply(Result.StatusData,function(x) dict.solveStatus[dict.solveStatus[,1] == x ,2])
  Result.status<-data.frame(code = Result.StatusData, desc = Result.StatusDesc)
  rownames(Result.status) <- dmu




  return(
    list(
      "rts" = Result.rts,
      "Orientation" = Result.orientation,
      "InputValues" =
        Result.Inputvalues,
      "OutputValues" = Result.Outputvalues,
      "Efficiency" = Result.eff,
      "Lambda" = Result.Lambda,
      "HCU_Input" = Hcu_input,
      "HCU_Output" = Hcu_output,
      "vx" = Result.vx,
      "uy" = Result.uy,
      "Free_Weights" = Result.Free_weights,
      "Model_Status" = Result.status
    )
  )
}

.Mult_OP <- function (x_t, y_t, rts = 'crs', weightRestriction) {
  Result.rts <- 'crs'

  if (rts == 'vrs') {
    Result.rts <- 'vrs'
  } else
  {
    Result.rts <- 'crs'
  }
  Result.orientation <- "Output"
  Result.Inputvalues <- x_t
  Result.Outputvalues <- y_t
  Result.eff_rec <- matrix(nrow = nrow(x_t))
  Result.Lambda <- matrix(nrow = nrow(x_t), ncol = nrow(x_t))
  Result.vx <- matrix(nrow = nrow(x_t), ncol = ncol(x_t))
  Result.uy <- matrix(nrow = nrow(y_t), ncol = ncol(y_t))
  Result.Free_weights <- matrix(nrow = nrow(x_t), ncol = 2)
  Result.StatusData <- vector(mode="logical", nrow(x_t))

  #free sign variable
  fs <- 0

  switch(rts,
         'crs' = fs <- 0,
         'vrs' = fs <- 2,
         fs <- 2)

  dmu <- row.names(x_t)

  # Assign names to Matrix for efficiency
  rownames(Result.eff_rec) <- dmu
  colnames(Result.eff_rec)[1] <- c("Eff")

  # Assign names to Matrix for efficiency
  rownames(Result.Lambda) <- dmu
  colnames(Result.Lambda) <- dmu

  # Assign names to Matrix for IP Weights
  rownames(Result.vx) <- dmu
  colnames(Result.vx) <- colnames(x_t)

  # Assign names to Matrix for OP Weights
  rownames(Result.uy) <- dmu
  colnames(Result.uy) <- colnames(y_t)

  # Assign names to Matrix for Free Weights
  rownames(Result.Free_weights) <- dmu
  colnames(Result.Free_weights) <- c("F1", "F2")


  for (d in 1:length(dmu)) {
    ip_v <- ncol(x_t)
    op_v <- ncol(y_t)

    # set objective function for DMU under consideration
    c <- 0
    for (i in 1:ncol(x_t)) {
      c[i] <- c(x_t[d, i])
    }

    #  for(i in 1:ncol(y_t)){
    # c[length(c) + 1] <- 0
    # }

    # If rts is vrs then add the free variable else omit
    if (rts == 'vrs') {
      # Adding the 1 and -1 for free variable
      if (length(c) != 0) {
        c[length(c) + 1] <- 1
        c[length(c) + 1] <- -1
      }
    }

    # indices for objective function
    ind_obj <- 0
    for (i in 1:ncol(x_t)) {
      ind_obj[i] <- i
    }

    # for(i in 1:ncol(y_t)){
    # ind_obj[ncol(x_t) + i ] <-nco1(x_t) + i
    # }

    # If rts is vrs then add the free variable indices else omit
    if (rts == 'vrs') {
      # indices for objective function fo r free variable
      if (length(ind_obj) != 0) {
        ind_obj[length(ind_obj) + 1] <- ip_v + op_v + 1
        ind_obj[length(ind_obj) + 1] <- ip_v + op_v + 2
      }
    }

    lps.model <- make.lp(0, ip_v + op_v + fs)
    set.objfn(lps.model, c, indices = ind_obj)


    # Output constraint for the DMU under consideration
    op_con_dmu <- 0
    for (i in 1:ncol(y_t)) {
      op_con_dmu[i] <- y_t[d, i]
    }

    #indices for DMU input constraints
    ind_op_con <- 0
    for (i in 1:ncol(y_t)) {
      ind_op_con[i] <- ip_v + i
    }
    add.constraint(lps.model, op_con_dmu, "=", rhs = 1, ind_op_con)

    # Multiply input values by -1 for building the constraints for all DMU's
    y_temp <- -1 * y_t

    # Construct constraints for all dmus
    for (i in 1:nrow(x_t)) {
      dmu_con <- 0
      for (j in 1:ncol(x_t)) {
        dmu_con[j] <- x_t[i, j]
      }
      for (k in 1:ncol(y_temp)) {
        dmu_con[ncol(x_t) + k] <- y_temp[i, k]
      }

      # if the rts is vrs then add free variable else omit
      if (rts == 'vrs') {
        if (length(dmu_con) != 0) {
          dmu_con[length(dmu_con) + 1] <- 1
          dmu_con[length(dmu_con) + 1] <- -1
        }
      }
      add.constraint(lps.model, dmu_con, ">=", rhs = 0)
    }

    #Weight Restrictions
    if (!missing(weightRestriction)) {
      colvalues <- c(colnames(x_t), colnames(y_t))

      for (i in 1:nrow(weightRestriction))
      {
        #First constraint
        add.constraint(lps.model,
                       c(weightRestriction[i, 1], -1),
                       ">=",
                       0,
                       indices = c(
                         which(colvalues == weightRestriction[i, 2]),
                         which(colvalues == weightRestriction[i, 3])
                       ))

        #Second constraint
        # Check if upper bound exists (Length ==4)
        if (length(weightRestriction) == 4) {
          #Check if upper bound is not Inf, NA, or NaN

          if (!is.nan(weightRestriction[i, 4]) &&
              !is.na(weightRestriction[i, 4]) &&
              !is.infinite(weightRestriction[i, 4])) {
            add.constraint(
              lps.model,
              c(1, -1 * weightRestriction[i, 4]),
              "<=",
              0,
              indices = c(
                which(colvalues == weightRestriction[i, 2]),
                which(colvalues == weightRestriction[i, 3])
              )
            )
          }
        }

      }
    }

    lp.control(lps.model, sense = 'min')


    Result.StatusData[d] <- solve(lps.model)

    Result.eff_rec[d] <- get.objective(lps.model)
    temp_lamdba <- get.dual.solution(lps.model)
    Result.Lambda[d, ] <- temp_lamdba[3:(length(dmu) + 2)]

    temp_weight <- get.variables(lps.model)
    #Store the IP & OP weights in respective Matrices
    Result.vx[d, 1:ncol(x_t)] <- temp_weight[1:ncol(x_t)]
    Result.uy[d, 1:ncol(y_t)] <-
      temp_weight[(ncol(x_t) + 1):(ncol(x_t) + ncol(y_t))]
    if (rts == 'vrs') {
      Result.Free_weights[d, 1:2] <-
        temp_weight[(ncol(x_t) + ncol(y_t) + 1):length(temp_weight)]
    }
  }
  Result.eff <- 1 / Result.eff_rec


  # Calculating the transpose of Lambda values for HCU calculation
  lambda_Transpose <- t(Result.Lambda)

  # Calculating HCU data for Input
  Hcu_input <- matrix(nrow = ncol(lambda_Transpose), ncol = ncol(x_t))
  colnames(Hcu_input) <- colnames(x_t)
  rownames(Hcu_input) <- rownames(lambda_Transpose)

  for (i in 1:ncol(lambda_Transpose)) {
    e <- x_t * lambda_Transpose[, i]
    for (j in 1:ncol(e)) {
      Hcu_input[i, j] <- sum(e[, j])
    }
  }

  # Calculating HCU data for Output
  Hcu_output <- matrix(nrow = ncol(lambda_Transpose), ncol = ncol(y_t))
  colnames(Hcu_output) <- colnames(y_t)
  rownames(Hcu_output) <- rownames(lambda_Transpose)

  for (i in 1:ncol(lambda_Transpose)) {
    e <- y_t * lambda_Transpose[, i]
    for (j in 1:ncol(e)) {
      Hcu_output[i, j] <- sum(e[, j])
    }
  }

  #Get the Status code and description from the LP model
  Result.StatusDesc<-sapply(Result.StatusData,function(x) dict.solveStatus[dict.solveStatus[,1] == x ,2])
  Result.status<-data.frame(code = Result.StatusData, desc = Result.StatusDesc)
  rownames(Result.status) <- dmu

  return(
    list(
      "rts" = Result.rts,
      "Orientation" = Result.orientation,
      "InputValues" =
        Result.Inputvalues,
      "OutputValues" = Result.Outputvalues,
      "Efficiency" = Result.eff,
      "Lambda" = Result.Lambda,
      "HCU_Input" = Hcu_input,
      "HCU_Output" = Hcu_output,
      "vx" = Result.vx,
      "uy" = Result.uy,
      "Free_Weights" = Result.Free_weights,
      "Model_Status" = Result.status
    )
  )
}

.Mal_Ben_Input<- function(x = x, y = y, rts ="crs",phase1_efficiency, phase = "mal", weightRestriction, include = TRUE){

  dmu <- row.names(x)

  Result.CEM <- matrix(nrow = nrow(x))
  rownames(Result.CEM) <- dmu
  colnames(Result.CEM)[1] <- c("Eff")

  Result.Lambda <- matrix(nrow = nrow(x), ncol = nrow(x))
  # Assign names to Matrix for efficiency
  rownames(Result.Lambda) <- dmu
  colnames(Result.Lambda) <- dmu

  Result.vx <-matrix(nrow = nrow(x), ncol = ncol(x))
  # Assign names to Matrix for IP Weights
  rownames(Result.vx) <- dmu
  colnames(Result.vx) <- colnames(x)

  Result.uy <-matrix(nrow = nrow(y), ncol = ncol(y))
  # Assign names to Matrix for OP Weights
  rownames(Result.uy) <- dmu
  colnames(Result.uy) <- colnames(y)

  Result.Free_weights <-matrix(nrow = nrow(x), ncol = 2)
  # Assign names to Matrix for Free Weights
  rownames(Result.Free_weights) <- dmu
  colnames(Result.Free_weights) <- c("F1","F2")

  Result.StatusData <- vector(mode="logical", nrow(x))

  sense <- ''
  equality <- ''

  # Set the Equality and equality based on Phase
  if(phase == 'mal'){
    sense <- 'min'
    equality <- '<='
  }else if (phase == 'ben'){
    sense <- 'max'
    equality <- '<='
  }

  #free sign variable
  fs<-0

  switch(rts,
         'crs' = fs <-0,
         'vrs' = fs <-2,
         fs<-2)

  #Second phase iteration over DMU efficiency scores
  #Setting input and output for each dmu
  for(d in 1 : length(dmu)){

    ip_v <- ncol(x)
    op_v <- ncol(y)

    # Set objective function for VIrtual DMU which is sum of the output variables for all DMUs
    c<-0
    for(i in 1:ncol(y)){

      if(include == TRUE){
        c[i]<-c(sum(y[,i]))
      }else if (include == FALSE) {
        c[i] <-c(sum(y[,i]) - y[d,i])
      }
    }


    # if vrs add the free variables else omit
    if(rts == 'vrs'){
      # Adding the 1 and -1 for free variable
      if(length(c) != 0){
        c[length(c) + 1] <- 1
        c[length(c) + 1] <- -1
      }
    }

    # indices for objective function
    ind_obj<-0
    for(i in 1:ncol(y)){
      ind_obj[i] <-c(ip_v + i)
    }

    # if vrs add the free variables indices else omit
    if(rts == 'vrs') {
      # indices for objective function for free variable
      if(length(ind_obj) !=0){
        ind_obj[length(ind_obj) + 1] <- ind_obj[length(ind_obj)] +1
        ind_obj[length(ind_obj) + 1] <- ind_obj[length(ind_obj)] +1
      }
    }

    lps.model <- make.lp(0, ip_v + op_v + fs)
    set.objfn(lps.model, c, indices = ind_obj)

    # Input constraint for the Virtual DMU which is sum of the Input variables for all DMUs
    ip_con_dmu<-0
    for(i in 1:ncol(x)){
      if(include == TRUE){
        ip_con_dmu[i] <- c(sum(x[,i]))
      }else if (include == FALSE) {
        ip_con_dmu[i] <- c(sum(x[,i]) - x[d,i])
      }
    }

    #indices for DMU input constraints
    ind_ip_con<-0
    for(i in 1:ncol(x)){
      ind_ip_con[i]<-i
    }

    add.constraint(lps.model, ip_con_dmu, "=", rhs = 1,ind_ip_con)

    # Multiply input values by -1 for building the constraints for all OMU's
    x_temp <- -1 *x

    # Construct constraints for all dmus excluding Rating DMU
    for(i in 1:nrow(x)){
      dmu_con<-0

      if(i != d){

        for(j in 1:ncol(x_temp)){
          dmu_con[j]<-x_temp[i,j]
        }
        for(k in 1:ncol(y)){
          dmu_con[ncol(x_temp) + k] <- y[i,k]
        }

        # If rts is vrs then add the free variable else omit
        if(rts == 'vrs') {
          #Adding free variable in the constraint
          if(length(dmu_con) !=0 ) {
            dmu_con[length(dmu_con) + 1] <- 1
            dmu_con[length(dmu_con) + 1] <- -1
          }
        }
        add.constraint(lps.model, dmu_con, equality, rhs = 0)
      }


      # Add constrains for the rating DMU
      if(i == d){
        dmu_con_Rating<-0
        for(j in 1:ncol(x_temp)){
          dmu_con_Rating[j]<- phase1_efficiency[d] * x_temp[d,j]
        }
        for(k in 1:ncol(y)){
          dmu_con_Rating[ncol(x_temp) + k] <- y[d,k]
        }

        # If rts is vrs then add the free variable else omit
        if(rts == 'vrs') {
          #Adding free variable in the constraint
          if(length(dmu_con_Rating) !=0 ) {
            dmu_con_Rating[length(dmu_con_Rating) + 1] <- 1
            dmu_con_Rating[length(dmu_con_Rating) + 1] <- -1
          }
        }
        add.constraint(lps.model, dmu_con_Rating, "=", rhs = 0)

      }
    }

    #Weight Restrictions
    if (!missing(weightRestriction)) {
      colvalues <- c(colnames(x), colnames(y))

      for (i in 1:nrow(weightRestriction))
      {
        #First constraint
        add.constraint(lps.model,
                       c(weightRestriction[i, 1], -1),
                       ">=",
                       0,
                       indices = c(
                         which(colvalues == weightRestriction[i, 2]),
                         which(colvalues == weightRestriction[i, 3])
                       ))
        #Second constraint
        # Check if upper bound exists (Length ==4)
        if (length(weightRestriction) == 4) {
          #Check if upper bound is not Inf, NA, or NaN

          if (!is.nan(weightRestriction[i, 4]) &&
              !is.na(weightRestriction[i, 4]) &&
              !is.infinite(weightRestriction[i, 4])) {
            add.constraint(
              lps.model,
              c(1, -1 * weightRestriction[i, 4]),
              "<=",
              0,
              indices = c(
                which(colvalues == weightRestriction[i, 2]),
                which(colvalues == weightRestriction[i, 3])
              )
            )
          }
        }
      }
    }

    lp.control(lps.model,sense=sense)

    Result.StatusData[d]<-solve(lps.model)


    Result.CEM[d]<-get.objective(lps.model)
    temp_lamdba <- get.dual.solution(lps.model)
    Result.Lambda[d,] <- temp_lamdba[3:(length(dmu)+2)]


    temp_weight <- get.variables(lps.model)

    ##Store the IP & OP weights in respective Matrices
    Result.vx[d,1:ncol(x)] <- temp_weight[1:ncol(x)]
    Result.uy[d,1:ncol(y)] <- temp_weight[(ncol(x)+1):(ncol(x) + ncol(y))]
    if(rts == 'vrs'){
      Result.Free_weights[d,1:2] <- temp_weight[(ncol(x) + ncol(y) + 1) : length(temp_weight)]
    }
  }

##Second Phase Cross Efficiency Calculation
  #Create the Cross Efficiency Matrix
  CrossEvaluation_Matrix <-
    matrix(nrow = length(dmu), ncol = length(dmu))
  rownames(CrossEvaluation_Matrix) <- dmu
  colnames(CrossEvaluation_Matrix) <- dmu

  #create Cross efficiency Average
  ce_ave <- matrix(nrow = 1, ncol = length(dmu))
  colnames(ce_ave) <- dmu
  rownames(ce_ave) <- "Average"

  #create Cross efficiency Max
  ceva_max <- matrix(nrow = 1, ncol = length(dmu))
  colnames(ceva_max) <- dmu
  rownames(ceva_max) <- "Max"

  #create Cross efficiency Min
  ceva_min <- matrix(nrow = 1, ncol = length(dmu))
  colnames(ceva_min) <- dmu
  rownames(ceva_min) <- "Min"


  # Input oriented CR efficiency
  for (i in 1:length(dmu)) {
    for (j in 1:length(dmu)) {
      #Calculate for IP values of DMU i * Ip weight of other DMUs J
      cr_ip <- x[i, ] * Result.vx[j, ]

      #Calculate for OP values of DMU i * OP weight of other DMUs J
      cr_op <- y[i, ] * Result.uy[j, ]

      #Calculate and store the values in Cross efficiency Matrix
      CrossEvaluation_Matrix[j, i] <- sum(cr_op) / sum(cr_ip)

    }
  }

  #Calculate cross efficiency Average
  for (i in 1:length(dmu)) {
    ce_ave[, i] <- mean(CrossEvaluation_Matrix[, i])
  }

  #Calculate cross efficiency Max
  for (i in 1:length(dmu)) {
    ceva_max[, i] <- max(CrossEvaluation_Matrix[, i])
  }

  #Calculate cross efficiency Max
  for (i in 1:length(dmu)) {
    ceva_min[, i] <- min(CrossEvaluation_Matrix[, i])
  }


  #Get the Status code and description from the LP model
  Result.StatusDesc<-sapply(Result.StatusData,function(x) dict.solveStatus[dict.solveStatus[,1] == x ,2])
  Result.status<-data.frame(code = Result.StatusData, desc = Result.StatusDesc)
  rownames(Result.status) <- dmu

  return(list("Phase2_Efficiency" = Result.CEM, "Phase2_Lambda" = Result.Lambda, "Phase2_vx" = Result.vx, "Phase2_uy" = Result.uy, "Phase2_Free_weights" = Result.Free_weights,"ceva_matrix" = CrossEvaluation_Matrix,
              "ce_ave" = ce_ave, "ceva_max" = ceva_max, "ceva_min" = ceva_min, "Phase2_Model_Status" = Result.status))
}

.Mal_Ben_Output<- function(x = x, y = y, rts ="crs",phase1_efficiency, phase = "mal", weightRestriction, include = TRUE){

  dmu <- row.names(x)

  Result.CEM <- matrix(nrow = nrow(x))
  rownames(Result.CEM) <- dmu
  colnames(Result.CEM)[1] <- c("Eff")

  Result.Lambda <- matrix(nrow = nrow(x), ncol = nrow(x))
  # Assign names to Matrix for efficiency
  rownames(Result.Lambda) <- dmu
  colnames(Result.Lambda) <- dmu

  Result.vx <-matrix(nrow = nrow(x), ncol = ncol(x))
  # Assign names to Matrix for IP Weights
  rownames(Result.vx) <- dmu
  colnames(Result.vx) <- colnames(x)

  Result.uy <-matrix(nrow = nrow(y), ncol = ncol(y))
  # Assign names to Matrix for OP Weights
  rownames(Result.uy) <- dmu
  colnames(Result.uy) <- colnames(y)

  Result.Free_weights <-matrix(nrow = nrow(x), ncol = 2)
  # Assign names to Matrix for Free Weights
  rownames(Result.Free_weights) <- dmu
  colnames(Result.Free_weights) <- c("F1","F2")

  Result.StatusData <- vector(mode="logical", nrow(x))

  sense <- ''
  equality <- ''

  # Set the Equality and equality based on Phase
  if(phase == 'mal'){
    sense <- 'min'
    equality <- '>='
  }else if (phase == 'ben'){
    sense <- 'max'
    equality <- '>='
  }

  #free sign variable
  fs<-0

  switch(rts,
         'crs' = fs <-0,
         'vrs' = fs <-2,
         fs<-2)

  #Second phase iteration over DMU efficiency scores
  #Setting input and output for each dmu
  for(d in 1:length(dmu)){

    ip_v <- ncol(x)
    op_v <- ncol(y)

    # set objective function for Virtual DMU which is sum of the output variables for all DMUs
    c<-0
    for(i in 1:ncol(x)){
      if(include == TRUE){
        c[i]<- c(sum(x[,i]))
      }else if (include == FALSE){
        c[i]<- c(sum(x[,i]) - x[d,i])
      }
    }


    # If rts is vrs then add the free variable else omit
    if(rts == 'vrs'){
      # Adding the 1 and -1 for free variable
      if(length(c) !=0){
        c[length(c) + 1] <- 1
        c[length(c) + 1] <- -1
      }
    }

    # indices for objective function
    ind_obj<-0
    for(i in 1:ncol(x)){
      ind_obj[i] <- i
    }


    # If rts is vrs then add the free variable indices else omit
    if(rts == 'vrs'){
      # indices for objective function fo r free variable
      if(length(ind_obj) != 0){
        ind_obj[length(ind_obj) + 1] <- ip_v + op_v + 1
        ind_obj[length(ind_obj) + 1] <- ip_v + op_v + 2
      }
    }

    lps.model <- make.lp(0, ip_v + op_v + fs)
    set.objfn(lps.model, c, indices = ind_obj)


    # Output constraint for the Virtual DMU which is sum of the Input variables for all DMUs
    op_con_dmu<-0
    for(i in 1:ncol(y)){
      if(include == TRUE){
        op_con_dmu[i] <- c(sum(y[,i]))
      }else if (include == FALSE) {
        op_con_dmu[i] <- c(sum(y[,i]) - y[d,i])
      }
    }

    #indices for DMU output constraints
    ind_op_con<-0
    for(i in 1:ncol(y)){
      ind_op_con[i]<- ip_v + i
    }
    add.constraint(lps.model, op_con_dmu, "=", rhs = 1,ind_op_con)

    # Multiply output values by -1 for building the constraints for all DMU's
    y_temp <- -1 * y

    # Construct constraints for all dmus excluding Rating DMU
    for(i in 1:nrow(x)){
      dmu_con<-0

      if(i != d){

        for(j in 1:ncol(x)){
          dmu_con[j] <-x[i,j]
        }
        for(k in 1:ncol(y_temp)){
          dmu_con[ncol(x) + k] <- y_temp[i,k]
        }

        # if the rts is vrs then add free variable else omit
        if(rts == 'vrs'){
          if(length(dmu_con) != 0 ) {
            dmu_con[length(dmu_con) + 1 ] <- 1
            dmu_con[length(dmu_con) + 1] <- -1
          }
        }
        add.constraint(lps.model, dmu_con, equality, rhs = 0)
      }

      # Add constrains for the rating DMU

      if(i == d){
        dmu_con_Rating<-0

        for(j in 1:ncol(x)){
          dmu_con_Rating[j] <-x[d,j]
        }

        for(k in 1:ncol(y_temp)){
          dmu_con_Rating[ncol(x) + k] <- phase1_efficiency[d] * y_temp[d,k]
        }

        # if the rts is vrs then add free variable else omit
        if(rts == 'vrs'){
          #Adding free variable in the constraint
          if(length(dmu_con_Rating) != 0 ) {
            dmu_con_Rating[length(dmu_con_Rating) + 1 ] <- 1
            dmu_con_Rating[length(dmu_con_Rating) + 1] <- -1
          }
        }
        add.constraint(lps.model, dmu_con_Rating, "=", rhs = 0)

      }
    }

    #Weight Restrictions
    if (!missing(weightRestriction)) {
      colvalues <- c(colnames(x), colnames(y))

      for (i in 1:nrow(weightRestriction))
      {
        #First constraint
        add.constraint(lps.model,
                       c(weightRestriction[i, 1], -1),
                       ">=",
                       0,
                       indices = c(
                         which(colvalues == weightRestriction[i, 2]),
                         which(colvalues == weightRestriction[i, 3])
                       ))

        #Second constraint
        # Check if upper bound exists (Length ==4)
        if (length(weightRestriction) == 4) {
          #Check if upper bound is not Inf, NA, or NaN

          if (!is.nan(weightRestriction[i, 4]) &&
              !is.na(weightRestriction[i, 4]) &&
              !is.infinite(weightRestriction[i, 4])) {
            add.constraint(
              lps.model,
              c(1, -1 * weightRestriction[i, 4]),
              "<=",
              0,
              indices = c(
                which(colvalues == weightRestriction[i, 2]),
                which(colvalues == weightRestriction[i, 3])
              )
            )
          }
        }

      }
    }

    lp.control(lps.model,sense=sense)
    Result.StatusData[d]<-solve(lps.model)


    Result.CEM[d]<- get.objective(lps.model)
    temp_lamdba <- get.dual.solution(lps.model)
    Result.Lambda[d,] <- temp_lamdba[3:(length(dmu)+2)]

    temp_weight <- get.variables(lps.model)

    #Store the IP & OP weights in respective Matrices
    Result.vx[d,1:ncol(x)] <- temp_weight[1:ncol(x)]
    Result.uy[d,1:ncol(y)] <- temp_weight[(ncol(x)+1):(ncol(x) + ncol(y))]
    if(rts == 'vrs'){
      Result.Free_weights[d,1:2] <- temp_weight[(ncol(x) + ncol(y) + 1) : length(temp_weight)]
    }
  }

  ##Second Phase Cross Efficiency Calculation
  #Create the Cross Efficiency Matrix
  CrossEvaluation_Matrix <-
    matrix(nrow = length(dmu), ncol = length(dmu))
  rownames(CrossEvaluation_Matrix) <- dmu
  colnames(CrossEvaluation_Matrix) <- dmu

  #create Cross efficiency Average
  ce_ave <- matrix(nrow = 1, ncol = length(dmu))
  colnames(ce_ave) <- dmu
  rownames(ce_ave) <- "Average"

  #create Cross efficiency Max
  ceva_max <- matrix(nrow = 1, ncol = length(dmu))
  colnames(ceva_max) <- dmu
  rownames(ceva_max) <- "Max"

  #create Cross efficiency Min
  ceva_min <- matrix(nrow = 1, ncol = length(dmu))
  colnames(ceva_min) <- dmu
  rownames(ceva_min) <- "Min"

  # Output oriented CR efficiency
  for (i in 1:length(dmu)) {
    for (j in 1:length(dmu)) {
      #Calculate for IP values of DMU i * Ip weight of other DMUs J
      cr_ip <- x[i, ] * Result.vx[j, ]

      #Calculate for OP values of DMU i * OP weight of other DMUs J
      cr_op <- y[i, ] * Result.uy[j, ]

      #Calculate and store the values in Cross efficiency Matrix
      CrossEvaluation_Matrix[j, i] <- 1 / (sum(cr_ip) / sum(cr_op))

    }
  }


  #Calculate cross efficiency Average
  for (i in 1:length(dmu)) {
    ce_ave[, i] <- mean(CrossEvaluation_Matrix[, i])
  }

  #Calculate cross efficiency Max
  for (i in 1:length(dmu)) {
    ceva_max[, i] <- max(CrossEvaluation_Matrix[, i])
  }

  #Calculate cross efficiency Max
  for (i in 1:length(dmu)) {
    ceva_min[, i] <- min(CrossEvaluation_Matrix[, i])
  }

  #Get the Status code and description from the LP model
  Result.StatusDesc<-sapply(Result.StatusData,function(x) dict.solveStatus[dict.solveStatus[,1] == x ,2])
  Result.status<-data.frame(code = Result.StatusData, desc = Result.StatusDesc)
  rownames(Result.status) <- dmu

  return(list("Phase2_Efficiency" = 1/Result.CEM, "Phase2_Lambda" = Result.Lambda, "Phase2_vx" = Result.vx, "Phase2_uy" = Result.uy, "Phase2_Free_weights" = Result.Free_weights, "ceva_matrix" = CrossEvaluation_Matrix,
              "ce_ave" = ce_ave, "ceva_max" = ceva_max, "ceva_min" = ceva_min, "Phase2_Model_Status" = Result.status))
}
