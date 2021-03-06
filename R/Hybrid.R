#' Calculate optimal number of clusters.
#'
#' This function calculates the optimal number of clusters for a given dataset.
#'
#' @param data data matrix
#' @param Treatment Vector of length p, indicating replicates of different treatment groups. For example, \emph{Treatment} = c(1,1,2,2,3,3) indicates 3 treatment groups, each with 2 replicates.
#'
#' @return A positive integer indicating the optimal number of clusters
#'
#' @export
#' @importFrom stats cor optimize pchisq quantile rmultinom
#'
#' @examples
#' Hybrid(sample_data, rep(c(1,2,3), each = 10))
#'
Hybrid = function(data, Treatment){
  K0 = floor(sqrt(nrow(data)))
  i0 = c(match(unique(Treatment), Treatment), length(Treatment) + 1)
  mydata = find_norm(data, Treatment)
  s = mydata$Normalizer
  dis = dis_tau(mydata)

  starting = initial(mydata, K0, dis)
  q0 = starting$q0
  mu0 = starting$mu0

  fn = hp_cluster(mydata, q0, mu0, method = 'EM')

  final = fn$final
  mu = fn$mu
  q = fn$q

  steps = rep(0, K0-1)
  r0 = rep(0, K0-1)

  nn = c()

  try = function(i, j, final){
    dataij = data[(final %in% c(i,j)), ]
    alphaij0 = alpha[(final %in% c(i,j))]
    a = g(dataij, s, alphaij0, Treatment, i0)
    alphaij = a$alpha
    muij = a$mu
    qij = a$q
    l = a$l

    return(list(alpha = alphaij, mu = muij, q = qij, l = l))
  }


  for(ind in 1:K0){
    #print(ind)
    if(ind == 1){
      # initial lglk's
      l0 = c()
      alpha = rep(0, nrow(data))

      for(i in 1:length(unique(final))){
        datai = data[final == i, ]
        if(is.vector(datai)) datai = matrix(datai, nrow = 1)
        else datai = as.matrix(datai, ncol = length(Treatment))
        alphai = (Est_alpha(s, datai, mu, Treatment))[, i]
        a = g(datai, s, alphai, Treatment, i0)
        alpha[final == i] = a$alpha

        l0 = c(l0, a$l)
      }

      l = matrix(1e40, nrow = length(l0), ncol = length(l0))
      for(i in 2:length(l0)){
        for(j in 1:(i-1)){
          l[i,j] = l0[i] + l0[j] - try(i, j, final)$l
        }
      }

    }

    #print(final)
    #print(l)
    min_l = min(l)
    merg = as.vector(which(l == min_l, arr.ind = TRUE))
    #print(merg)
    k1 = which(merg[1] == levels(final))
    k2 = which(merg[2] == levels(final))
    r = nrow(data[(final %in% merg), ])

    nn = c(nn, merg[2])

    l[merg[2], ] = 1e40
    l[, merg[2]] = 1e40

    steps[length(steps) - ind + 1] = min_l
    r0[length(r0) - ind + 1] = r

    mu[k1,] = try(merg[1], merg[2], final)$mu
    mu = mu[-k2, ]
    final[final %in% merg] = merg[1]
    final = as.factor(as.numeric(as.character(final)))

    if (length(levels(final)) == 1) break

    l0[merg[2]] = 0
    i = merg[1]
    datai = data[final == i, ]
    if(is.vector(datai)){datai = matrix(datai, nrow = 1)} else{datai = as.matrix(datai, ncol = length(Treatment))}
    alphai = (Est_alpha(s, datai, mu, Treatment))[, which(i == as.numeric(levels(final)))]
    a = g(datai, s, alphai, Treatment, i0)
    l0[merg[1]] = a$l

    for(j in 1:(merg[1]-1)){
      if(sum(j == nn) == 0){
        l[merg[1], j] = l0[j] + l0[merg[1]] - try(j, merg[1], final)$l
      }
    }

    if(i != K0){
      for(j in (merg[1]+1):nrow(l)){
        if(sum(j == nn) == 0){l[j, merg[1]] = l0[j] + l0[merg[1]] - try(j, merg[1], final)$l}
      }
    }

  }

  h = c(steps, r0)
  I = length(unique(Treatment))

  p_v = c()
  for(i in 1:(length(h)/2)){
    temp = 1 - pchisq(h[i], df = h[i+length(h)/2] + 2*I-1)
    p_v = c(p_v, temp)
  }

  if(sum(p_v > 0.01) > 0) k = min(which(p_v > 0.01))
  else k = length(p_v) + 1

  return(k)
}




