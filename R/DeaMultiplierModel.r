#******************************************************************************
#
# Copyright Aurobindh Kalathil Puthanpura, ETA, PSU, 2015, 2016
# Use granted under BSD license terms
#
# R DeaMultiplierModel Package
#
# NOTE: See utility.R for description of naming nomenclature, use of .l & .b suffix
#
# Multiplier DEA wraper function
#
#******************************************************************************
#

DeaMultiplierModel<-function (x = x, y = y, rts = "crs", orientation = "input"){


  rts         <- .checkoption(rts,          "rts",          options.rts.l)
  orientation <- .checkoption(orientation,  "orientation",  options.orientation.l)


  # Check that x & y are legal inputs & convert to standard values
  x <- .checkData(x, "x")
  y <- .checkData(y, "y")

  .checkDataGood(x, y)

  if (nrow(x) != nrow(y))
    stop("Number of DMU's in inputs != number of DMU's in outputs", call. = FALSE)

  # Check the orientation and rts to decide the internal function to choose

  if(.loadLibrary()!=TRUE)
      stop("The required package lpSolveAPI could not be installed. Please install and try again", call. = FALSE)

  if(orientation == "input" && rts == "crs"){
    results <- .Mult_lP(x, y, "crs")
    return(results)
  }else if(orientation == "input" && rts == "vrs"){
    results <- .Mult_lP(x, y, "vrs")
    return(results)
  }else if(orientation == "output" && rts == "crs"){
    results <- .Mult_OP(x, y, "crs")
    return(results)
  }else if(orientation == "output" && rts == "vrs"){
    results <- .Mult_OP(x, y, "vrs")
    return(results)
  }
}

CrossEfficiency <-function (x = x, y = y, rts = "crs", orientation = "input"){

#Pass the input values to DeaMultiplierModel function and obtain the results
result <- DeaMultiplierModel(x,y,rts,orientation)

# Get the row names from the InputValues
dmu <-rownames(result$InputValues)
length(dmu)

#Create the Cross Efficiency Matrix
CrossEfficiency_Matrix <- matrix(nrow = length(dmu), ncol = length(dmu))
rownames(CrossEfficiency_Matrix) <- dmu
colnames(CrossEfficiency_Matrix) <- dmu

#create Cross efficiency Average
ce_ave <- matrix(nrow = 1, ncol = length(dmu))
colnames(ce_ave) <- dmu
rownames(ce_ave) <- "Average"

#create Cross efficiency Max
ce_max <- matrix(nrow = 1, ncol = length(dmu))
colnames(ce_max) <- dmu
rownames(ce_max) <- "Max"

#create Cross efficiency Min
ce_min <- matrix(nrow = 1, ncol = length(dmu))
colnames(ce_min) <- dmu
rownames(ce_min) <- "Min"

if(result$Orientation == 'Input'){
  # Input oriented CR efficiency
  for(i in 1:length(dmu)){
    for(j in 1:length(dmu)){
      #Calculate for IP values of DMU i * Ip weight of other DMUs J
      cr_ip <- result$InputValues[i,] * result$IP_Weights[j,]

      #Calculate for OP values of DMU i * OP weight of other DMUs J
      cr_op <- result$OutputValues[i,] * result$OP_Weights[j,]

      #Calculate and store the values in Cross efficiency Matrix
      CrossEfficiency_Matrix[j,i] <- sum(cr_op) /sum(cr_ip)

    }
  }
}else if(result$Orientation == 'Output'){
  # Input oriented CR efficiency
  for(i in 1:length(dmu)){
    for(j in 1:length(dmu)){
      #Calculate for IP values of DMU i * Ip weight of other DMUs J
      cr_ip <- result$InputValues[i,] * result$IP_Weights[j,]

      #Calculate for OP values of DMU i * OP weight of other DMUs J
      cr_op <- result$OutputValues[i,] * result$OP_Weights[j,]

      #Calculate and store the values in Cross efficiency Matrix
      CrossEfficiency_Matrix[j,i] <-sum(cr_ip)/sum(cr_op)

    }
  }
}

#Calculate cross efficiency Average
for(i in 1:length(dmu)){
  ce_ave[,i]<- mean(CrossEfficiency_Matrix[,i])
}

#Calculate cross efficiency Max
for(i in 1:length(dmu)){
  ce_max[,i]<- max(CrossEfficiency_Matrix[,i])
}

#Calculate cross efficiency Max
for(i in 1:length(dmu)){
  ce_min[,i]<- min(CrossEfficiency_Matrix[,i])
}

return(list("ce_matrix" = CrossEfficiency_Matrix, "ce_ave" = ce_ave, "ce_max" = ce_max, "ce_min" = ce_min))

}


#Mal_Ben<-function(x = x, y = y, rts ="crs", orientation = "input", phase = "Mal", include = TRUE){

#}
