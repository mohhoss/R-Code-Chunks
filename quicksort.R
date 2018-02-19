quicksort <- function(arr, p = 1, r = length(arr)) {
     
     if (!is.vector(arr)) {
          arr <- as.vector(arr)
     }
     
     arr <- quick_sort(arr , p, r)
     return(arr)
}

quick_sort <- function(arr, p ,r) {
     if (p < r) {
          output <- partition(arr, p , r)
          q <- output$Q
          arr <- output$A
          arr <- quick_sort(arr , p , r =q-1)
          arr <- quick_sort(arr , p = q+1 , r= r)
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

