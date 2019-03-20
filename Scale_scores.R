# Function to rescale risk scores to new range
# dat = matrix of calculated risk scores
# newMin = desired minimum score 
# newMax = desired maximum score 

new.range <- function (dat, newMin, newMax) {
  result <- sapply(dat, function (x) {
    ((newMax - newMin) * (x - min(dat))) / (max(dat) - min(dat)) + newMin 
  })
  return(matrix(result, nrow=nrow(dat)))
}
