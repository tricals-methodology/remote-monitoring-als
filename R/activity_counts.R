# Algorithm Activity Counts based on Neishabouri et al. (2022)
# Quantification of acceleration as activity counts in ActiGraph wearable.
# Sci Rep. 2022 Jul 13;12(1):11958.
# doi: 10.1038/s41598-022-16003-x. PMID: 35831446

# Code adapted from https://github.com/actigraph/agcounts
# Currently only works for data sampled at 30 Hz

# Define 7th order IIR filter
library (gsignal)

# Input data coefficients (b)
INPUT_COEFFICIENTS <- c(
  -0.009341062898525, -0.025470289659360, -0.004235264826105,
  0.044152415456420, 0.036493718347760, -0.011893961934740,
  -0.022917390623150, -0.006788163862310, 0.000000000000000
)


# Output data coefficients (a)
OUTPUT_COEFFICIENTS <- c(
  1.00000000000000000000, -3.63367395910957000000, 5.03689812757486000000,
  -3.09612247819666000000, 0.50620507633883000000, 0.32421701566682000000,
  -0.15685485875559000000, 0.01949130205890000000, 0.00000000000000000000
)


bpf_filter <- function(raw_data) {
  # Input: Raw accelerometer data 
  # Band-pass filter data (data has to have a 30 Hz sampling frequency; if not, downsample first)
  # data is rounded to 3 decimal places before input into BPF.
  raw_data <- round(raw_data, 3)

  # Compute initial condition for filter
  zi <- c(gsignal::filter_zi(INPUT_COEFFICIENTS, OUTPUT_COEFFICIENTS))

  # Compute filtered data
  bpf_data <- gsignal::filter(
    filt = INPUT_COEFFICIENTS,
    a = OUTPUT_COEFFICIENTS,
    x = raw_data,
    zi = zi * raw_data[1]
  )[[1]]

  # Clean up memory
  rm(raw_data)
  gc()

  # Scale data in order to replicate the range of the AM7164.
  bpf_data <- ((3.0 / 4096.0) / (2.6 / 256.0) * 237.5) * bpf_data

  # Return filtered data
  return(bpf_data)
}


trim_data <- function(bpf_data, lfe_select) {
  # Low Frequency Extension (LFE)
  if (lfe_select) {
    min_count <- 1
    max_count <- 128

    # Rectify the rescaled signal
    trim_data <- abs(bpf_data)
    trim_data[(trim_data < min_count) & (trim_data >= 4)] <- 0
    trim_data[trim_data > max_count] <- max_count
    mask <- (trim_data < 4) & (trim_data >= min_count)
    trim_data[mask] <- abs(trim_data[mask]) - 1
    trim_data <- floor(trim_data)
    rm(mask)
  } else {
    min_count <- 4
    max_count <- 128

    # Rectify the rescaled signal
    trim_data <- abs(bpf_data)
    trim_data[trim_data < min_count] <- 0
    trim_data[trim_data > max_count] <- max_count
    trim_data <- floor(trim_data)
  }

  return(trim_data)
}


# Further down-sample the signal and low-pass filter to 10 Hz by a non-overlapping moving average
resample_10hz <- function(trim_data) {
  # hackish downsample to 10 Hz
  downsample_10hz <- cumsum(trim_data)
  downsample_10hz[4:length(downsample_10hz)] <- downsample_10hz[4:length(downsample_10hz)] - downsample_10hz[1:(length(downsample_10hz) - 3)]
  downsample_10hz <- floor(downsample_10hz[seq(from = 3, to = length(downsample_10hz), by = 3)] / 3)

  return(downsample_10hz)
}


# Finally, find the counts by summing the down-sampled signal within the predefined epoch length for each axis:
sum_counts <- function(downsample_10hz, epoch_seconds) {
  # Accumulator for epoch
  block_size <- epoch_seconds * 10
  epoch_counts <- cumsum(downsample_10hz)

  epoch_counts[(block_size + 1):length(epoch_counts)] <- epoch_counts[(block_size + 1):length(epoch_counts)] - epoch_counts[1:(length(epoch_counts) - block_size)]

  extract_indices <- (block_size - 1) %% block_size + 1 + seq(0, length(epoch_counts) - block_size, by = block_size)
  epoch_counts <- epoch_counts[extract_indices]

  epoch_counts <- floor(epoch_counts)

  return(epoch_counts)
}








