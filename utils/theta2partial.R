theta2partial <- function(theta) {
  if(length(dim(theta))==2){
    d <- solve(sqrt(diag(diag(theta)+1e-12)))
    partial <- -d %*% theta %*% d
    diag(partial)<-1
    return(partial)
  }else if (length(dim(theta))==3){
    partial <- theta
    for (k in 1:dim(theta)[1]){
      d <- solve(sqrt(diag(diag(theta[k,,]))))
      partial[k,,] <- -d %*% theta[k,,] %*% d
      diag(partial[k,,])<-1
    }
    return(partial)
  }
  
}