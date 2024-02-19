# Algorithm classification non-wear times based on Choi et al. (2011)
# Validation of accelerometer wear and nonwear time classification algorithm. 
# Med Sci Sports Exerc. 2011 Feb;43(2):357-64. 
# doi: 10.1249/MSS.0b013e3181ed61a3. PMID: 20581716

# With optimized hyperparameters based on Syed et al. (2022)
# Evaluating the performance of raw and epoch non-wear algorithms using multiple accelerometers and electrocardiogram recordings. 
# Sci Rep. 2020 Apr 3;10(1):5866. doi: 10.1038/s41598-020-62821-2. PMID: 32246080

# Code adapted from https://github.com/shaheen-syed/ActiGraph-ActiWave-Analysis

library(PhysicalActivity)

## Description
## 1-min time intervals with consecutive zero counts for at least 90-min time window (window 1), 
# allowing a short time intervals with nonzero counts lasting up to 2 min (allowance interval) 
# if no counts are detected during both the 30 min (window 2) of upstream and downstream from that interval; 
# any nonzero counts except the allowed short interval are considered as wearing
choi_2011_calculate_non_wear_time <- function(data,
                                              activity_threshold = 0, # Minimum spike threshold (counts)
                                              min_period_len = 210, # Minimum interval (mins)
                                              spike_tolerance = 1, # Allowed #spikes
                                              min_window_len = 20, # Up/Down window size (mins)
                                              window_spike_tolerance = 0, # Allowed #spikes in window
                                              epoch_seconds = 10, # Epoch length (s)
                                              print_output = FALSE) {
  
  # Input: aggregated data in Activity Counts
  # Collapse data to 60s epoch length
  if (epoch_seconds != 60) {
    aggr <- PhysicalActivity::dataCollapser(data, TS = "DATE", by = 60, col = "VM")
  } else {
    aggr <- data
  }

  # Check if aggr contains at least min_period_len of aggr
  if (nrow(aggr) < min_period_len) {
    stop(paste("Epoch aggr contains", nrow(aggr), "samples, which is less than the", min_period_len, "minimum required samples"))
  }

  # Create non wear vector as numeric vector with ones. now we only need to add the zeros which are the non-wear time segments
  non_wear_vector <- rep(1, nrow(aggr))

  ##. Variables used to keep track of non wear periods
  # indicator for resetting and starting over
  reset <- FALSE
  # indicator for stopping the non-wear period
  stopped <- FALSE
  # indicator for starting to count the non-wear period
  start <- FALSE
  # second window validation
  window_2_invalid <- FALSE
  # starting minute for the non-wear period
  strt_nw <- 0
  # ending minute for the non-wear period
  end_nw <- 0
  # counter for the number of minutes with intensity between 1 and 100
  cnt_non_zero <- 0
  # keep track of non wear sequences
  ranges <- data.frame()

  # loop over the aggr
  for (n in 1:nrow(aggr)) {
    # get the value
    nval <- aggr[n, "VM"]

    # reset counters if reset or stopped
    if (reset == TRUE | stopped == TRUE) {
      strt_nw <- 0
      end_nw <- 0
      start <- FALSE
      reset <- FALSE
      stopped <- FALSE
      window_2_invalid <- FALSE
      cnt_non_zero <- 0
    }

    # the non-wear period starts with a zero count
    if (nval == 0 & start == FALSE) {
      # assign the starting minute of non-wear
      strt_nw <- n
      # set start boolean to true so we know that we started the period
      start <- TRUE
    }

    # only do something when the non-wear period has started
    if (start == TRUE) {
      # keep track of the number of minutes with intensity that is not a 'zero' count
      if (nval > activity_threshold) {
        # increase the spike counter
        cnt_non_zero <- cnt_non_zero + 1
      }

      # when there is a non-zero count, check the upstream and downstream window for counts
      # only when the upstream and downstream window have zero counts, then it is a valid non wear sequence
      if (nval > 0) {
        if (n == nrow(aggr)) {
          upstream <- 0
        } else {
          # check upstream window if there are counts, note that we skip the count right after the spike, since we allow for 2 minutes of spikes
          upstream <- aggr[(n + spike_tolerance):ifelse(
            n + min_window_len + (spike_tolerance - 1) > nrow(aggr),
            nrow(aggr) - (spike_tolerance - 1),
            n + min_window_len + (spike_tolerance - 1)
          ), "VM"]
        }

        # check if upstream has non zero counts, if so, then the window is invalid
        if (sum(upstream > 0) > window_spike_tolerance) {
          window_2_invalid <- TRUE
        }

        # check downstream window if there are counts, again, we skip the count right before since we allow for 2 minutes of spikes
        downstream <- aggr[(ifelse((n - min_window_len) > 0, n - min_window_len, 1)):(n - 1), "VM"]

        # check if downstream has non zero counts, if so, then the window is invalid
        if (sum(downstream > 0) > window_spike_tolerance) {
          window_2_invalid <- TRUE
        }

        # if the second window is invalid, we need to reset the sequence for the next run
        if (window_2_invalid) {
          reset <- TRUE
        }
      }


      if (nval <= activity_threshold) {
        cnt_non_zero <- 0
      }

      # the sequence ends when there are 3 consecutive spikes, or an invalid second window (upstream or downstream), or the last value of the sequence
      if (cnt_non_zero == (spike_tolerance + 1) || window_2_invalid == TRUE || n == nrow(aggr)) {
        # define the end of the period
        end_nw <- n

        # check if the sequence is sufficient in length
        if (nrow(aggr[strt_nw:end_nw, ]) < min_period_len) {
          # length is not sufficient, so reset values in next run
          reset <- TRUE
        } else {
          # length of sequence is sufficient, set stopped to TRUE so we save the sequence start and end later on
          stopped <- TRUE
        }
      }

      # if stopped is TRUE, the sequence stopped and is valid to include in the ranges
      if (stopped == T) {
        # add ranges start and end non wear time
        ranges <- rbind(ranges, c(strt_nw, end_nw))
      }
    }
  }

  # ignore when no non wear times were flagged
  if (nrow(ranges) > 0) {
    # set the non wear vector according to start and end
    for (r in 1:nrow(ranges)) {
      # set the non wear vector according to start and end
      non_wear_vector[ranges[r, 1]:(ranges[r, 2] - 1)] <- 0

      # if set to True, then print output to console/log
      if (print_output == TRUE) {
        print(paste(
          "No. no wear:", r, "start aggregated index:", ranges[r, 1], "end aggregated index:",
          ranges[r, 2], "length:", ranges[r, 2] - ranges[r, 1]
        ))
        print(paste(
          "start time:", aggr[ranges[r, 1], "DATE"], "end time:", aggr[ranges[r, 2], "DATE"],
          "duration in hrs:", round(as.numeric(aggr[ranges[r, 2], "DATE"] - aggr[ranges[r, 1], "DATE"]), digits = 2)
        ))
      }
    }
  }

  # translate back to epoch length
  non_wear_vector <- rep(non_wear_vector, each = 60 / epoch_seconds)

  data$WEAR <- non_wear_vector[1:nrow(data)]

  return(data)
}
