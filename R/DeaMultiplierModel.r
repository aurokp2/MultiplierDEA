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
        "ce_ave" = t(ce_ave),
        "ceva_max" = t(ceva_max),
        "ceva_min" = t(ceva_min),
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
