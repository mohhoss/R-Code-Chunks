quicksort <- function(arr, increasing = TRUE, randomized = TRUE) {
     
     # Initialization
     if (!is.vector(arr)) {
          arr <- as.vector(arr)
     }
     p = 1
     r = length(arr)
     
     # Call to the recursive function
     arr <- quick_sort(arr , p, r, rnd = randomized)
     
     # implementing the increasing/decreasing arrangement
     if (!increasing) {
          arr <- rev(arr)
     }
     
     return(arr)
}

# Main recursive function
quick_sort <- function(arr, p ,r, rnd) {
     if (p < r) {
          
          # implementing the randomized version of the algorithm
          if (rnd) {
               i <- sample(p:r, 1)
               temp <- arr[r]
               arr[r] <- arr[i]
               arr[i] <- temp
               rm(temp)
          }
          
          output <- partition(arr, p , r)
          q <- output$Q
          arr <- output$A
          arr <- quick_sort(arr , p , r =q-1, rnd)
          arr <- quick_sort(arr , p = q+1 , r= r, rnd)
     }
     return(arr)
}
partition <- function(arr, p, r) {

     x = arr[r]
     i = p -1
     
     for (j in p:(r-1)) {
          
          if (arr[j] <= x) {
               i = i+1
               temp <- arr[i]
               arr[i] <- arr[j]
               arr [j] <- temp
               rm(temp)
          }
     }
     i = i+1
     temp <- arr[i]
     arr[i] <- arr[r]
     arr[r] <- temp
     rm(temp)
     
     output <- list(Q = i, A =arr)
     return(output)
}
