#**************************************************************************
#
# Copyright Aurobindh Kalathil Puthanpura, ETA, PSU, 2015
# Use granted under BSD license terms
#
# R DEA M u l t ip l ie r Model Package
#
# Author: AUrobindh Kalathil Puthanpura
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
deaMM.version = "0.1.0"
deaMM.id = paste0("deaMM version ", deaMM.version)

# Load ipsolverAPI library

.loadLibrary <-function(){

  if("lpSolveAPI" %in% rownames(installed.packages()))
  {
    library("lpSolveAPI")
  }else
  {
    install.packages("lpSolveAPI")
    library("lpSolveAPI")
  }

  if("lpSolveAPI" %in% rownames(installed.packages()))
  {
    return(TRUE)
  }else
  {
    return(FALSE)
  }

}

#
# Options Strings
#

options.orientation.l <- c ("input", "output")
options.rts.l <- c ("vrs" , "crs")

# Function: .checkDate
# Check for valid types
# Check add names to unnamed variables

.checkData <- function(x, name){
# Is value a dataframe, array or vector?
if ( !is.data.frame(x) && !is.array (x ) && (length(x) < 2))
stop(name," is not a matrix, array (1 or 2 dimensions) or data.frame", call. =
FALSE)
if(length(dim(x)) > 2)
  stop(name," is greater then 2 dimensions", call. = FALSE)

# If data.frame - convert to array
if ( is.data.frame(x))
x <- data.matrix(x)

# If vector convert to a 1 x N array
if( is.vector (x ) )
X <- array(x, c ( length(x ), 1))

# Check that a l l numeric
if(!is.numeric(x))
stop(name," must be numeric", call. = FALSE)

# Add DMU names i f empty
if (length(dim(x)) >= 1 && is.null(rownames(x)))
rownames(x) <- paste0("DMU", l:nrow(x))

# Add Col names i f empty
if (length(dim(x)) >= 2 && is.null(colnames(x)))

colnames(x) <- paste0(name, l:ncol(x))
return(x)
}


# Function .CheckDataGood
# Check input & output data and warn about problems for dea
.checkDataGood <- function(x, y){

status <- TRUE

# Check for any zero's
# check for na' s *
# Check for no positive values on input & output
if (any(x == 0, na.rm=TRUE)){
  cat("Warning, data has DMU's with inputs that are zero, th is may cause numerical
  problems\n")
  status <- FALSE
}

if (any(y == 0, na.rm=TRUE)){
cat("Warning, data has DMU's with outputs that are zero, this may cause
numerical problems\n")
status <- FALSE
}

for (k in (1 :nrow (x))) {
  if (all(x[k,] <= 0, na.rm=TRUE) & all(y[k,] <= 0, na.rm=TRUE)){
    cat("warning, dmu # " ,k , "has no positive inputs or outputs, this may cause
    numerical problems\n")
    status <- FALSE
  }
}
return(status)
}

# Function .checkoption()
# Check that value in legal l i s t options
# If on legal list, return lowercase short form option
#
.checkoption <- function(value, name, options.l=c(true, false)){
  if (length(value) != 1){
       stop("illegal vector for ",name, " option must be single value", call. = false)
  }
             # If options are character
             if (is.character(options.l[1])){
             if (is.character(value)){
             tmp.value <- tolower(value)
             # if (tmp.value %in% options.1)
             i <- charmatch(tmp.value, options.l, nomatch=-1)
             if ( i > 0)
              return(options.l[i])
             }
             } else if (is.logical(Options.l[1])){
             # Logical options
             if (is.logical(value)) {
             return(value)
             }
             } else {
# Numeric options
  if(is.numeric(value) && is.finite(value) && (value >= 0)){
    return(value)
  }
  options.l <- c("numeric >= 0")
}
stop("Illegal value=", value, " for ", name, "; legal values are: ",
  paste(options.l, collapse = ", " ) ,
  call. = FALSE)
}

.Mult_lP <- function (x_t, y_t, rts = 'crs'){
Result.rts <- 'crs'
if(rts == 'vrs'){
Result.rts <- 'vrs'
} else
{
Result.rts <- 'crs'
}

Result.orientation <-"Input"
Result.Inputvalues <- x_t
Result.Outputvalues <-y_t
Result.eff <- matrix(nrow = nrow(x_t))
Result.Lambda <- matrix(nrow = nrow(x_t), ncol = nrow(x_t))
Result.IP_Weights <-matrix(nrow = nrow(x_t), ncol = ncol(x_t))
Result.OP_Weights <-matrix(nrow = nrow(y_t), ncol = ncol(y_t))
Result.Free_weights <-matrix(nrow = nrow(x_t), ncol = 2)


#free sign variable
fs<-0

switch(rts,
        'crs' = fs <-0,
        'vrs' = fs <-2,
         fs<-2)

dmu <- row.names(x_t)

# Assign names to Matrix for efficiency
rownames(Result.eff) <- dmu
colnames(Result.eff)[1] <-c("Eff")

# Assign names to Matrix for efficiency
rownames(Result.Lambda) <- dmu
colnames(Result.Lambda) <- dmu

# Assign names to Matrix for IP Weights
rownames(Result.IP_Weights) <- dmu
colnames(Result.IP_Weights) <- colnames(x_t)

# Assign names to Matrix for OP Weights
rownames(Result.OP_Weights) <- dmu
colnames(Result.OP_Weights) <- colnames(y_t)

# Assign names to Matrix for Free Weights
rownames(Result.Free_weights) <- dmu
colnames(Result.Free_weights) <- c("F1","F2")



for(d in 1 : length(dmu)){

ip_v <- ncol(x_t)
op_v <- ncol(y_t)

# Set objective function for DMU under consideration
c<-0
for(i in 1:ncol(y_t)){
  c[i]<-c(y_t[d,i])
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
for(i in 1:ncol(y_t)){
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

# Input constraint for the DMU under consideration
ip_con_dmu<-0
for(i in 1:ncol(x_t)){
  ip_con_dmu[i] <-x_t[d,i]
}

#indices for DMU input constraints
ind_ip_con<-0
for(i in 1:ncol(x_t)){
    ind_ip_con[i]<-i
}

add.constraint(lps.model, ip_con_dmu, "=", rhs = 1,ind_ip_con)

# Multiply input values by -1 for building the constraints for all OMU's
x_temp <- -1 *x_t

# Construct constraints for all dmus
for(i in 1:nrow(x_t)){
    dmu_con<-0
    for(j in 1:ncol(x_temp)){
        dmu_con[j]<-x_temp[i,j]
    }
     for(k in 1:ncol(y_t)){
          dmu_con[ncol(x_temp) + k] <- y_t[i,k]
     }

    # If rts is vrs then add the free variable else omit
    if(rts == 'vrs') {
      #Adding free variable in the constraint
      if(length(dmu_con) !=0 ) {
         dmu_con[length(dmu_con) + 1] <- 1
         dmu_con[length(dmu_con) + 1] <- -1
      }

    }
    add.constraint(lps.model, dmu_con, "<=", rhs = 0)
}
lp.control(lps.model,sense='max')

solve(lps.model)

Result.eff[d]<-get.objective(lps.model)
temp_lamdba <- get.dual.solution(lps.model)
Result.Lambda[d,] <- temp_lamdba[3:(length(dmu)+2)]


temp_weight <- get.variables(lps.model)

#Store the IP & OP weights in respective Matrices
Result.IP_Weights[d,1:ncol(x_t)] <- temp_weight[1:ncol(x_t)]
Result.OP_Weights[d,1:ncol(y_t)] <- temp_weight[(ncol(x_t)+1):(ncol(x_t) + ncol(y_t))]
if(rts == 'vrs'){
Result.Free_weights[d,1:2] <- temp_weight[(ncol(x_t) + ncol(y_t) + 1) : length(temp_weight)]
}
}

# Calculating the transpose of Lambda values for HCU calculation
lambda_Transpose <- t(Result.Lambda)

# Calculating HCU data for Input
Hcu_input <-matrix(nrow = ncol(lambda_Transpose), ncol = ncol(x_t))
colnames(Hcu_input) <- colnames(x_t)
rownames(Hcu_input) <- rownames(lambda_Transpose)

for(i in 1:ncol(lambda_Transpose)) {
  e <- x_t * lambda_Transpose[,i]
  for(j in 1:ncol(e)){
    Hcu_input[i,j] <- sum(e[,j])
  }
}

# Calculating HCU data for Output
Hcu_output <-matrix(nrow = ncol(lambda_Transpose), ncol = ncol(y_t))
colnames(Hcu_output) <- colnames(y_t)
rownames(Hcu_output) <- rownames(lambda_Transpose)

for(i in 1:ncol(lambda_Transpose)) {
  e <- y_t * lambda_Transpose[,i]
  for(j in 1:ncol(e)){
    Hcu_output[i,j] <- sum(e[,j])
  }
}

return(list("rts" =Result.rts, "Orientation" = Result.orientation, "InputValues" =
Result.Inputvalues, "OutputValues" = Result.Outputvalues, "Efficiency" = Result.eff, "Lambda" = Result.Lambda,
"HCU_Input" = Hcu_input, "HCU_Output" = Hcu_output, "IP_Weights" = Result.IP_Weights, "OP_Weights" = Result.OP_Weights, "Free_Weights" = Result.Free_weights))
}

.Mult_OP <- function (x_t, y_t, rts = 'crs'){
Result.rts <- 'crs'

if(rts == 'vrs'){
Result.rts <- 'vrs'
} else
{
Result.rts <- 'crs'
}
Result.orientation <-"Output"
Result.Inputvalues <- x_t
Result.Outputvalues <-y_t
Result.eff_rec <- matrix(nrow = nrow(x_t))
Result.Lambda <- matrix(nrow = nrow(x_t), ncol = nrow(x_t))
Result.IP_Weights <-matrix(nrow = nrow(x_t), ncol = ncol(x_t))
Result.OP_Weights <-matrix(nrow = nrow(y_t), ncol = ncol(y_t))
Result.Free_weights <-matrix(nrow = nrow(x_t), ncol = 2)

#free sign variable
fs<-0

switch(rts,
      'crs' = fs <-0,
      'vrs' = fs <-2,
      fs<-2)

dmu <- row.names(x_t)

# Assign names to Matrix for efficiency
rownames(Result.eff_rec) <- dmu
colnames(Result.eff_rec)[1] <-c("Eff")

# Assign names to Matrix for efficiency
rownames(Result.Lambda) <- dmu
colnames(Result.Lambda) <- dmu

# Assign names to Matrix for IP Weights
rownames(Result.IP_Weights) <- dmu
colnames(Result.IP_Weights) <- colnames(x_t)

# Assign names to Matrix for OP Weights
rownames(Result.OP_Weights) <- dmu
colnames(Result.OP_Weights) <- colnames(y_t)

# Assign names to Matrix for Free Weights
rownames(Result.Free_weights) <- dmu
colnames(Result.Free_weights) <- c("F1","F2")


for(d in 1:length(dmu)){

  ip_v <- ncol(x_t)
  op_v <- ncol(y_t)

# set objective function for DMU under consideration
c<-0
for(i in 1:ncol(x_t)){
  c[i]<- c(x_t[d,i])
}

#  for(i in 1:ncol(y_t)){
# c[length(c) + 1] <- 0
# }

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
for(i in 1:ncol(x_t)){
  ind_obj[i] <- i
}

# for(i in 1:ncol(y_t)){
# ind_obj[ncol(x_t) + i ] <-nco1(x_t) + i
# }

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


# Output constraint for the DMU under consideration
op_con_dmu<-0
for(i in 1:ncol(y_t)){
  op_con_dmu[i] <- y_t[d,i]
}

#indices for DMU input constraints
ind_op_con<-0
for(i in 1:ncol(y_t)){
  ind_op_con[i]<- ip_v + i
}
add.constraint(lps.model, op_con_dmu, "=", rhs = 1,ind_op_con)

# Multiply input values by -1 for building the constraints for all DMU's
y_temp <- -1 * y_t

# Construct constraints for all dmus
for(i in 1:nrow(x_t)){
    dmu_con<-0
    for(j in 1:ncol(x_t)){
      dmu_con[j] <-x_t[i,j]
}
for(k in 1:ncol(y_temp)){
  dmu_con[ncol(x_t) + k] <- y_temp[i,k]
}

# if the rts is vrs then add free variable else omit
  if(rts == 'vrs'){
    if(length(dmu_con) != 0 ) {
      dmu_con[length(dmu_con) + 1 ] <- 1
        dmu_con[length(dmu_con) + 1] <- -1
    }
  }
   add.constraint(lps.model, dmu_con, ">=", rhs = 0)
}

lp.control(lps.model,sense='min')

solve(lps.model)

Result.eff_rec[d]<- get.objective(lps.model)
temp_lamdba <- get.dual.solution(lps.model)
Result.Lambda[d,] <- temp_lamdba[3:(length(dmu)+2)]

temp_weight <- get.variables(lps.model)
#Store the IP & OP weights in respective Matrices
Result.IP_Weights[d,1:ncol(x_t)] <- temp_weight[1:ncol(x_t)]
Result.OP_Weights[d,1:ncol(y_t)] <- temp_weight[(ncol(x_t)+1):(ncol(x_t) + ncol(y_t))]
if(rts == 'vrs'){
  Result.Free_weights[d,1:2] <- temp_weight[(ncol(x_t) + ncol(y_t) + 1) : length(temp_weight)]
}
}
Result.eff <- 1/Result.eff_rec


# Calculating the transpose of Lambda values for HCU calculation
lambda_Transpose <- t(Result.Lambda)

# Calculating HCU data for Input
Hcu_input <-matrix(nrow = ncol(lambda_Transpose), ncol = ncol(x_t))
colnames(Hcu_input) <- colnames(x_t)
rownames(Hcu_input) <- rownames(lambda_Transpose)

for(i in 1:ncol(lambda_Transpose)) {
  e <- x_t * lambda_Transpose[,i]
  for(j in 1:ncol(e)){
    Hcu_input[i,j] <- sum(e[,j])
  }
}

# Calculating HCU data for Output
Hcu_output <-matrix(nrow = ncol(lambda_Transpose), ncol = ncol(y_t))
colnames(Hcu_output) <- colnames(y_t)
rownames(Hcu_output) <- rownames(lambda_Transpose)

for(i in 1:ncol(lambda_Transpose)) {
  e <- y_t * lambda_Transpose[,i]
  for(j in 1:ncol(e)){
    Hcu_output[i,j] <- sum(e[,j])
  }
}

return(list("rts" = Result.rts, "Orientation" = Result.orientation, "InputValues" =
           Result.Inputvalues, "OutputValues" = Result.Outputvalues, "Efficiency" = Result.eff, "Lambda" = Result.Lambda,
           "HCU_Input" = Hcu_input, "HCU_Output" = Hcu_output, "IP_Weights" = Result.IP_Weights, "OP_Weights" = Result.OP_Weights, "Free_Weights" = Result.Free_weights))
}

