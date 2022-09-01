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

DeaMultiplierModel <-
  function (x = x,
            y = y,
            rts = "crs",
            orientation = "input",
            weightRestriction) {
    rts         <-
      .checkoption(rts,          "rts",          options.rts.l)
    orientation <-
      .checkoption(orientation,  "orientation",  options.orientation.l)


    # Check that x & y are legal inputs & convert to standard values
    x <- .checkData(x, "x")
    y <- .checkData(y, "y")

    .checkDataGood(x, y)

    if (!missing(weightRestriction)) {
      .checkWeightRestriction(x, y, weightRestriction)
    }

    if (nrow(x) != nrow(y))
      stop("Number of DMU's in inputs != number of DMU's in outputs", call. = FALSE)

    # if (.loadLibrary() != TRUE)
    #   stop(
    #     "The required package lpSolveAPI could not be installed. Please install and try again",
    #     call. = FALSE
    #   )

    # Check the orientation and rts to decide the internal function to choose

    if (orientation == "input" && rts == "crs") {
      results <- .Mult_lP(x, y, "crs", weightRestriction)
      return(results)
    } else if (orientation == "input" && rts == "vrs") {
      results <- .Mult_lP(x, y, "vrs", weightRestriction)
      return(results)
    } else if (orientation == "output" && rts == "crs") {
      results <- .Mult_OP(x, y, "crs", weightRestriction)
      return(results)
    } else if (orientation == "output" && rts == "vrs") {
      results <- .Mult_OP(x, y, "vrs", weightRestriction)
      return(results)
    }
  }

CrossEfficiency <-
  function (x = x,
            y = y,
            rts = "crs",
            orientation = "input",
            weightRestriction) {
    #Pass the input values to DeaMultiplierModel function and obtain the results
    result <- DeaMultiplierModel(x, y, rts, orientation, weightRestriction)

    # Get the row names from the InputValues
    dmu <- rownames(result$InputValues)
    length(dmu)

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

    if (result$Orientation == 'Input') {
      # Input oriented CR efficiency
      for (i in 1:length(dmu)) {
        for (j in 1:length(dmu)) {
          #Calculate for IP values of DMU i * Ip weight of other DMUs J
          cr_ip <- result$InputValues[i, ] * result$vx[j, ]

          #Calculate for OP values of DMU i * OP weight of other DMUs J
          cr_op <- result$OutputValues[i, ] * result$uy[j, ]

          #Calculate and store the values in Cross efficiency Matrix
          CrossEvaluation_Matrix[j, i] <- sum(cr_op) / sum(cr_ip)

        }
      }
    } else if (result$Orientation == 'Output') {
      # Output oriented CR efficiency
      for (i in 1:length(dmu)) {
        for (j in 1:length(dmu)) {
          #Calculate for IP values of DMU i * Ip weight of other DMUs J
          cr_ip <- result$InputValues[i, ] * result$vx[j, ]

          #Calculate for OP values of DMU i * OP weight of other DMUs J
          cr_op <- result$OutputValues[i, ] * result$uy[j, ]

          #Calculate and store the values in Cross efficiency Matrix
          CrossEvaluation_Matrix[j, i] <- 1 / (sum(cr_ip) / sum(cr_op))

        }
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

    return(
      list(
        "ceva_matrix" = CrossEvaluation_Matrix,
        "ce_ave" = ce_ave,
        "ceva_max" = ceva_max,
        "ceva_min" = ceva_min,
        "vx" = result$vx,
        "uy" = result$uy,
        "Model_Status" = result$Model_Status
      )
    )

  }


Mal_Ben <-
  function(x = x,
           y = y,
           rts = "crs",
           orientation = "input",
           phase = "mal",
           weightRestriction,
           include = TRUE) {
    rts         <-
      .checkoption(rts,          "rts",          options.rts.l)
    orientation <-
      .checkoption(orientation,  "orientation",  options.orientation.l)
    phase         <-
      .checkoption(phase,          "phase",          options.phase.l)


    #Pass the input values to DeaMultiplierModel function and obtain the results
    result <- DeaMultiplierModel(x, y, rts, orientation, weightRestriction)

    #Pass result from first phase to second phase function
    if (orientation == "input" && rts == "crs" && phase == "mal") {
      results <-
        .Mal_Ben_Input(x, y, "crs", result$Efficiency, "mal", weightRestriction, include)
      return(
        list(
          "rts" = result$rts,
          "Orientation" = result$Orientation,
          "InputValues" = result$InputValues,
          "OutputValues" = result$OutputValues,
          "Phase1_Efficiency" = result$Efficiency,
          "Phase1_Lambda" = result$Lambda,
          "Phase1_vx" = result$vx,
          "Phase1_uy" = result$uy,
          "Phase1_Free_Weights" = result$Free_Weights,
          "Phase1_Model_Status" = result$Model_Status,
          "Phase2_Efficiency" = results$Phase2_Efficiency,
          "Phase2_Lambda" = results$Phase2_Lambda,
          "Phase2_vx" = results$Phase2_vx,
          "Phase2_uy" = results$Phase2_uy,
          "Phase2_Free_weights" = results$Phase2_Free_weights,
          "Phase2_Model_Status" = results$Phase2_Model_Status,
          "ceva_matrix" = results$ceva_matrix,
          "ce_ave" = results$ce_ave,
          "ceva_max" = results$ceva_max,
          "ceva_min" = results$ceva_min
        )
      )

    } else if (orientation == "input" && rts == "crs" && phase == "ben") {
      results <-
        .Mal_Ben_Input(x, y, "crs", result$Efficiency, "ben", weightRestriction, include)
      return(
        list(
          "rts" = result$rts,
          "Orientation" = result$Orientation,
          "InputValues" = result$InputValues,
          "OutputValues" = result$OutputValues,
          "Phase1_Efficiency" = result$Efficiency,
          "Phase1_Lambda" = result$Lambda,
          "Phase1_vx" = result$vx,
          "Phase1_uy" = result$uy,
          "Phase1_Free_Weights" = result$Free_Weights,
          "Phase1_Model_Status" = result$Model_Status,
          "Phase2_Efficiency" = results$Phase2_Efficiency,
          "Phase2_Lambda" = results$Phase2_Lambda,
          "Phase2_vx" = results$Phase2_vx,
          "Phase2_uy" = results$Phase2_uy,
          "Phase2_Free_weights" = results$Phase2_Free_weights,
          "Phase2_Model_Status" = results$Phase2_Model_Status,
          "ceva_matrix" = results$ceva_matrix,
          "ce_ave" = results$ce_ave,
          "ceva_max" = results$ceva_max,
          "ceva_min" = results$ceva_min
        )
      )

    } else if (orientation == "input" && rts == "vrs" && phase == "mal") {
      results <-
        .Mal_Ben_Input(x, y, "vrs", result$Efficiency, "mal", weightRestriction, include)
      return(
        list(
          "rts" = result$rts,
          "Orientation" = result$Orientation,
          "InputValues" = result$InputValues,
          "OutputValues" = result$OutputValues,
          "Phase1_Efficiency" = result$Efficiency,
          "Phase1_Lambda" = result$Lambda,
          "Phase1_vx" = result$vx,
          "Phase1_uy" = result$uy,
          "Phase1_Free_Weights" = result$Free_Weights,
          "Phase1_Model_Status" = result$Model_Status,
          "Phase2_Efficiency" = results$Phase2_Efficiency,
          "Phase2_Lambda" = results$Phase2_Lambda,
          "Phase2_vx" = results$Phase2_vx,
          "Phase2_uy" = results$Phase2_uy,
          "Phase2_Free_weights" = results$Phase2_Free_weights,
          "Phase2_Model_Status" = results$Phase2_Model_Status,
          "ceva_matrix" = results$ceva_matrix,
          "ce_ave" = results$ce_ave,
          "ceva_max" = results$ceva_max,
          "ceva_min" = results$ceva_min
        )
      )

    } else if (orientation == "input" && rts == "vrs" && phase == "ben") {
      results <-
        .Mal_Ben_Input(x, y, "vrs", result$Efficiency, "ben", weightRestriction, include)
      return(
        list(
          "rts" = result$rts,
          "Orientation" = result$Orientation,
          "InputValues" = result$InputValues,
          "OutputValues" = result$OutputValues,
          "Phase1_Efficiency" = result$Efficiency,
          "Phase1_Lambda" = result$Lambda,
          "Phase1_vx" = result$vx,
          "Phase1_uy" = result$uy,
          "Phase1_Free_Weights" = result$Free_Weights,
          "Phase1_Model_Status" = result$Model_Status,
          "Phase2_Efficiency" = results$Phase2_Efficiency,
          "Phase2_Lambda" = results$Phase2_Lambda,
          "Phase2_vx" = results$Phase2_vx,
          "Phase2_uy" = results$Phase2_uy,
          "Phase2_Free_weights" = results$Phase2_Free_weights,
          "Phase2_Model_Status" = results$Phase2_Model_Status,
          "ceva_matrix" = results$ceva_matrix,
          "ce_ave" = results$ce_ave,
          "ceva_max" = results$ceva_max,
          "ceva_min" = results$ceva_min
        )
      )
      #---------- Malevolent and Benevolent for Output oriented
    } else if (orientation == "output" &&
               rts == "crs" && phase == "mal") {
      results <-
        .Mal_Ben_Output(x, y, "crs", result$Efficiency, "mal", weightRestriction, include)
      return(
        list(
          "rts" = result$rts,
          "Orientation" = result$Orientation,
          "InputValues" = result$InputValues,
          "OutputValues" = result$OutputValues,
          "Phase1_Efficiency" = result$Efficiency,
          "Phase1_Lambda" = result$Lambda,
          "Phase1_vx" = result$vx,
          "Phase1_uy" = result$uy,
          "Phase1_Free_Weights" = result$Free_Weights,
          "Phase1_Model_Status" = result$Model_Status,
          "Phase2_Efficiency" = results$Phase2_Efficiency,
          "Phase2_Lambda" = results$Phase2_Lambda,
          "Phase2_vx" = results$Phase2_vx,
          "Phase2_uy" = results$Phase2_uy,
          "Phase2_Free_weights" = results$Phase2_Free_weights,
          "Phase2_Model_Status" = results$Phase2_Model_Status,
          "ceva_matrix" = results$ceva_matrix,
          "ce_ave" = results$ce_ave,
          "ceva_max" = results$ceva_max,
          "ceva_min" = results$ceva_min
        )
      )

    } else if (orientation == "output" &&
               rts == "crs" && phase == "ben") {
      results <-
        .Mal_Ben_Output(x, y, "crs", result$Efficiency, "ben", weightRestriction, include)
      return(
        list(
          "rts" = result$rts,
          "Orientation" = result$Orientation,
          "InputValues" = result$InputValues,
          "OutputValues" = result$OutputValues,
          "Phase1_Efficiency" = result$Efficiency,
          "Phase1_Lambda" = result$Lambda,
          "Phase1_vx" = result$vx,
          "Phase1_uy" = result$uy,
          "Phase1_Free_Weights" = result$Free_Weights,
          "Phase1_Model_Status" = result$Model_Status,
          "Phase2_Efficiency" = results$Phase2_Efficiency,
          "Phase2_Lambda" = results$Phase2_Lambda,
          "Phase2_vx" = results$Phase2_vx,
          "Phase2_uy" = results$Phase2_uy,
          "Phase2_Free_weights" = results$Phase2_Free_weights,
          "Phase2_Model_Status" = results$Phase2_Model_Status,
          "ceva_matrix" = results$ceva_matrix,
          "ce_ave" = results$ce_ave,
          "ceva_max" = results$ceva_max,
          "ceva_min" = results$ceva_min
        )
      )

    } else if (orientation == "output" &&
               rts == "vrs" && phase == "mal") {
      results <-
        .Mal_Ben_Output(x, y, "vrs", result$Efficiency, "mal", weightRestriction, include)
      return(
        list(
          "rts" = result$rts,
          "Orientation" = result$Orientation,
          "InputValues" = result$InputValues,
          "OutputValues" = result$OutputValues,
          "Phase1_Efficiency" = result$Efficiency,
          "Phase1_Lambda" = result$Lambda,
          "Phase1_vx" = result$vx,
          "Phase1_uy" = result$uy,
          "Phase1_Free_Weights" = result$Free_Weights,
          "Phase1_Model_Status" = result$Model_Status,
          "Phase2_Efficiency" = results$Phase2_Efficiency,
          "Phase2_Lambda" = results$Phase2_Lambda,
          "Phase2_vx" = results$Phase2_vx,
          "Phase2_uy" = results$Phase2_uy,
          "Phase2_Free_weights" = results$Phase2_Free_weights,
          "Phase2_Model_Status" = results$Phase2_Model_Status,
          "ceva_matrix" = results$ceva_matrix,
          "ce_ave" = results$ce_ave,
          "ceva_max" = results$ceva_max,
          "ceva_min" = results$ceva_min
        )
      )

    } else if (orientation == "output" &&
               rts == "vrs" && phase == "ben") {
      results <-
        .Mal_Ben_Output(x, y, "vrs", result$Efficiency, "ben", weightRestriction, include)
      return(
        list(
          "rts" = result$rts,
          "Orientation" = result$Orientation,
          "InputValues" = result$InputValues,
          "OutputValues" = result$OutputValues,
          "Phase1_Efficiency" = result$Efficiency,
          "Phase1_Lambda" = result$Lambda,
          "Phase1_vx" = result$vx,
          "Phase1_uy" = result$uy,
          "Phase1_Free_Weights" = result$Free_Weights,
          "Phase1_Model_Status" = result$Model_Status,
          "Phase2_Efficiency" = results$Phase2_Efficiency,
          "Phase2_Lambda" = results$Phase2_Lambda,
          "Phase2_vx" = results$Phase2_vx,
          "Phase2_uy" = results$Phase2_uy,
          "Phase2_Free_weights" = results$Phase2_Free_weights,
          "Phase2_Model_Status" = results$Phase2_Model_Status,
          "ceva_matrix" = results$ceva_matrix,
          "ce_ave" = results$ce_ave,
          "ceva_max" = results$ceva_max,
          "ceva_min" = results$ceva_min
        )
      )

    }

  }

SDEA <- function(x=x, y=y, orientation="input", rts ="crs", Cook=FALSE){
  # check ip, op, orientation and rts
  rts <- .checkoption(rts,          "rts",          options.rts.l)
  orientation <- .checkoption(orientation,  "orientation",  options.orientation.l)


  # Check that x & y are legal inputs & convert to standard values
  x <- .checkData(x, "x")
  y <- .checkData(y, "y")

  .checkDataGood(x, y)

  # Check the orientation and rts to decide the internal function to choose
  if (Cook==FALSE){
    results <- .sdea_internal(x=x, y=y,orientation=orientation, rts=rts)
    return(results)
  }else if(Cook==TRUE){
    results <- .sdea_cook_internal(x=x, y=y,orientation=orientation, rts=rts,cook=TRUE)
    return(results)
  }

}


MPI <-function(Dataset=Dataset, DMU_colName=DMU_colName, IP_colNames=IP_colNames, OP_ColNames=OP_ColNames, Period_ColName=Period_ColName, Periods=Periods, rts="crs", orientation="input", scale=FALSE){

  rts <- .checkoption(rts,          "rts",          options.rts.l)
  orientation <- .checkoption(orientation,  "orientation",  options.orientation.l)

  # Dataframe to store all
  results <- data.frame(stringsAsFactors = FALSE)

  # Compare for all the years
  for(j in 1:(length(Periods)-1)){
    # store temporary results
    results_temp <- data.frame(stringsAsFactors = FALSE)

    # Get data for two time periods
    data_t1 <- Dataset %>% filter(eval(parse(text = Period_ColName)) == Periods[j])
    data_t2 <- Dataset %>% filter(eval(parse(text = Period_ColName)) == Periods[j + 1])

    # Get Dmu names from time period t
    dmu_t1 <- data_t1[,DMU_colName]

    # IP columns for period t
    x_temp_t1 <- as.data.frame(data_t1[IP_colNames])
    row.names(x_temp_t1) <- dmu_t1

    # OP columns for period t
    y_temp_t1 <- as.data.frame(data_t1[OP_ColNames])
    row.names(y_temp_t1) <- dmu_t1


    # Get Dmu names from time period t + 1
    dmu_t2 <- data_t2[,DMU_colName]

    # IP columns for period t + 1
    x_temp_t2 <- as.data.frame(data_t2[IP_colNames])
    row.names(x_temp_t2) <- dmu_t2

    # OP columns for period t + 1
    y_temp_t2 <- as.data.frame(data_t2[OP_ColNames])
    row.names(y_temp_t2) <- dmu_t2

    # Assign temporary dataframe with Unique DMU names from t and t + 1
    results_temp <- c(as.character(dmu_t1), as.character(dmu_t2)) %>% as.data.frame()
    colnames(results_temp)[1] <- "DMU"
    results_temp <- results_temp %>% select("DMU") %>% unique()

    # Efficiency for DMU in t with reference to time period t
    # Determine crs et1t1
    if (rts == "crs" || scale){
      et1t1_temp<-.sdea_mpi_internal(x=x_temp_t1, y=y_temp_t1, orientation=orientation, rts="crs", Cook=FALSE, ref.x=x_temp_t1, ref.y=y_temp_t1)
      et1t1 <- cbind(rownames(et1t1_temp$Eff), et1t1_temp$Eff) %>% as.data.frame()
      colnames(et1t1)[1] <- "DMU"
      colnames(et1t1)[2] <- "et1t1.crs"

      results_temp <- left_join(results_temp,et1t1, by = c("DMU" = "DMU"))
    }

    # Determine vrs et1t1
    if (rts == "vrs"){
      et1t1_temp<-.sdea_cook_mpi_internal(x=x_temp_t1, y=y_temp_t1, orientation=orientation, rts="vrs", Cook=TRUE, ref.x=x_temp_t1, ref.y=y_temp_t1)
      et1t1 <- cbind(rownames(et1t1_temp$Eff), et1t1_temp$Eff) %>% as.data.frame()
      colnames(et1t1)[1] <- "DMU"
      colnames(et1t1)[2] <- "et1t1.vrs"

      results_temp <- left_join(results_temp,et1t1, by = c("DMU" = "DMU"))
    }

    # Efficiency for DMU in t + 1 with reference to time period t + 1
    # Determine crs et2t2
    if (rts == "crs" || scale){
      et2t2_temp<-.sdea_mpi_internal(x=x_temp_t2, y=y_temp_t2, orientation=orientation, rts="crs", Cook=FALSE, ref.x=x_temp_t2, ref.y=y_temp_t2)
      et2t2 <- cbind(rownames(et2t2_temp$Eff), et2t2_temp$Eff) %>% as.data.frame()
      colnames(et2t2)[1] <- "DMU"
      colnames(et2t2)[2] <- "et2t2.crs"

      results_temp <- left_join(results_temp,et2t2, by = c("DMU" = "DMU"))
    }

    # Determine vrs et2t2
    if (rts == "vrs"){
      et2t2_temp<-.sdea_cook_mpi_internal(x=x_temp_t2, y=y_temp_t2, orientation=orientation, rts="vrs", Cook=TRUE, ref.x=x_temp_t2, ref.y=y_temp_t2)
      et2t2 <- cbind(rownames(et2t2_temp$Eff), et2t2_temp$Eff) %>% as.data.frame()
      colnames(et2t2)[1] <- "DMU"
      colnames(et2t2)[2] <- "et2t2.vrs"

      results_temp <- left_join(results_temp,et2t2, by = c("DMU" = "DMU"))
    }

    # Efficiency for DMU in t with reference to time period t + 1
    # Determine crs et1t2
    if (rts == "crs" || scale){
      et1t2_temp<-.sdea_mpi_internal(x=x_temp_t2, y=y_temp_t2, orientation=orientation, rts="crs", Cook=FALSE, ref.x=x_temp_t1, ref.y=y_temp_t1)
      et1t2 <- cbind(rownames(et1t2_temp$Eff), et1t2_temp$Eff) %>% as.data.frame()
      colnames(et1t2)[1] <- "DMU"
      colnames(et1t2)[2] <- "et1t2.crs"

      results_temp <- left_join(results_temp,et1t2, by = c("DMU" = "DMU"))
    }

    # Determine vrs et1t2
    if (rts == "vrs"){
      et1t2_temp<-.sdea_cook_mpi_internal(x=x_temp_t2, y=y_temp_t2, orientation=orientation, rts="vrs", Cook=TRUE, ref.x=x_temp_t1, ref.y=y_temp_t1)
      et1t2 <- cbind(rownames(et1t2_temp$Eff), et1t2_temp$Eff) %>% as.data.frame()
      colnames(et1t2)[1] <- "DMU"
      colnames(et1t2)[2] <- "et1t2.vrs"

      results_temp <- left_join(results_temp,et1t2, by = c("DMU" = "DMU"))
    }

    # Efficiency for DMU in t + 1 with reference to time period t
    # Determine crs et2t1
    if (rts == "crs" || scale){
      et2t1_temp<-.sdea_mpi_internal(x=x_temp_t1, y=y_temp_t1, orientation=orientation, rts="crs", Cook=FALSE, ref.x=x_temp_t2, ref.y=y_temp_t2)
      et2t1 <- cbind(rownames(et2t1_temp$Eff), et2t1_temp$Eff) %>% as.data.frame()
      colnames(et2t1)[1] <- "DMU"
      colnames(et2t1)[2] <- "et2t1.crs"

      results_temp <- left_join(results_temp,et2t1, by = c("DMU" = "DMU"))
    }

    # Determine vrs et2t1
    if (rts == "vrs"){
      et2t1_temp<-.sdea_cook_mpi_internal(x=x_temp_t1, y=y_temp_t1, orientation=orientation, rts="vrs", Cook=TRUE, ref.x=x_temp_t2, ref.y=y_temp_t2)
      et2t1 <- cbind(rownames(et2t1_temp$Eff), et2t1_temp$Eff) %>% as.data.frame()
      colnames(et2t1)[1] <- "DMU"
      colnames(et2t1)[2] <- "et2t1.vrs"

      results_temp <- left_join(results_temp,et2t1, by = c("DMU" = "DMU"))
    }

    # Calculate (pure) technical efficiency change
    if(rts == "crs") {
      results_temp["tec"] <- as.numeric(as.character(results_temp$et2t2.crs))/as.numeric(as.character(results_temp$et1t1.crs))
    }else{
      results_temp["ptec"] <- as.numeric(as.character(results_temp$et2t2.vrs))/as.numeric(as.character(results_temp$et1t1.vrs))
    }

    # Calculate scale efficiency change if using vrs
    if(rts == "vrs" && scale){
      results_temp["sec1"] <- (as.numeric(as.character(results_temp$et1t2.crs))/as.numeric(as.character(results_temp$et1t2.vrs)))/(as.numeric(as.character(results_temp$et1t1.crs))/as.numeric(as.character(results_temp$et1t1.vrs)))
      results_temp["sec2"] <- (as.numeric(as.character(results_temp$et2t2.crs))/as.numeric(as.character(results_temp$et2t2.vrs)))/(as.numeric(as.character(results_temp$et2t1.crs))/as.numeric(as.character(results_temp$et2t1.vrs)))
      results_temp["sec"] <- (results_temp$sec1 * results_temp$sec2) ^ 0.5
    }else
    {
      results_temp["sec1"] <- NA
      results_temp["sec2"] <- NA
      results_temp["sec"] <- NA
    }

    # Calculate technology change (technology frontier shift)
    if(rts == "crs"){
      results_temp["tc1"] <- as.numeric(as.character(results_temp$et1t2.crs))/as.numeric(as.character(results_temp$et2t2.crs))
      results_temp["tc2"] <- as.numeric(as.character(results_temp$et1t1.crs))/as.numeric(as.character(results_temp$et2t1.crs))
      results_temp["tc"] <- (results_temp$tc1 * results_temp$tc2) ^ 0.5
    }else {
      results_temp["tc1"] <- as.numeric(as.character(results_temp$et1t2.vrs))/as.numeric(as.character(results_temp$et2t2.vrs))
      results_temp["tc2"] <- as.numeric(as.character(results_temp$et1t1.vrs))/as.numeric(as.character(results_temp$et2t1.vrs))
      results_temp["tc"] <- (results_temp$tc1 * results_temp$tc2) ^ 0.5
    }

    # Calculate Malmquist index
    if(rts == "crs"){
      results_temp["m.crs"] <- results_temp$tec * results_temp$tc
    }else
    {
      results_temp["m.vrs"] <- results_temp$ptec * results_temp$tc
    }

    # Create new colume for Year  (t - t + 1)
    results_temp["Year"] <- paste(Periods[j], Periods[j + 1], sep="-")
    results <- rbind(results, results_temp)

  }
  row.names(results) <- seq(1:nrow(results))

  return(results)
}
