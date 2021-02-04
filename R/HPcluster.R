#' Hurdle Poisson clustering
#'
#' This function gives the clustering result based on a hurdle Poisson model.
#'
#' @param data Data matrix with dimension N*P
#' @param Treatment Vector of length P. Indicating replicates of different treatment groups. For example, \emph{Treatment} = c(1,1,2,2,3,3) indicates 3 treatment groups, each with 2 replicates.
#' @param nK Positive integer. Number of clusters.
#' @param method Method for the algorithm. Can choose from \emph{"EM"} as Expectation Maximization or \emph{"SA"} as Simulated Annealing.
#' @param cool Real number between (0, 1). Cooling rate for the \emph{"SA"} algorithm. Uses 0.9 by default.
#' @param nstart Positive integer. Number of starts for the entire algorithm. Note that as \emph{nstart} increases the computational time also grows linearly. Uses 1 by default.
#'
#' @return
#' \describe{
#' \item{cluster}{Vector consists of integers from 1 to nK, indicating final clustering result. For evaluating the clustering result please check \link[aricode]{NMI} for \emph{normalized mutual information}.}
#' \item{prob}{N*nK matrix, the (i, j)th element representing the probability that observation i belongs to cluster j.}
#' \item{log_l}{Scaler, the hurdle Poisson log-likelihood of the final clustering result.}
#' }
#'
#' @importFrom stats cor optimize pchisq quantile rmultinom
#' @export
#'
#' @examples HPcluster(sample_data, rep(c(1,2,3), 10), nK = 7, method = 'EM', cool = 0.9, nstart = 1)
HPcluster = function(data, Treatment, nK, method = c('EM', 'SA'), cool = 0.9, nstart = 1){
  mydata = find_norm(data, Treatment)

  l0 = c()
  f0 = list()
  Z0 = list()

  dis = dis_tau(mydata)

  for(tr in 1:nstart){
    starting = initial(mydata, nK, dis)
    q0 = starting$q0
    mu0 = starting$mu0
    fn = hp_cluster(mydata, q0, mu0, method = method, cool = cool)
    l0 = c(l0, fn$lglk)
    f0[[tr]] = fn$final
    Z0[[tr]] = fn$Z
  }

  f1 = f0[[which.max(l0)]]
  Z1 = Z0[[which.max(l0)]]
  l1 = max(l0)
  return(list(prob = Z1, cluster = f1, log_l = l1))
}
