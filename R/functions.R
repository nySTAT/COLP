require(combinat)
require(MASS)

COLP_exhaustive = function(y,x){
  nlev_y = nlevels(y)
  P =  combinat::permn(nlev_y)
  nP = length(P)
  xCy = rep(0,nP)
  for (l in 1:nP){
    yy=y
    levels(yy)=as.character(P[[l]])
    yy=as.factor(as.numeric(as.character(yy)))
    if (nlevels(yy)>2){
      xCy[l]=stats::logLik(MASS::polr(yy~x,method="logistic"))
    }else{
      xCy[l]=stats::logLik(stats::glm(yy~x,family=stats::binomial))
    }
  }
  i = which.max(xCy)
  return(list(M=xCy[i],P=P[[i]]))
}

COLP_greedy = function(y,x,P=NULL){
  nlev_y = nlevels(y)
  if (is.null(P)){
    P =  sample(nlev_y)
    #P = 1:nlev_y
  }
  Ps = P
  yy=y
  levels(yy)=as.character(P)
  nly=nlevels(yy)
  yy=as.factor(as.numeric(as.character(yy)))
  if (nly>2){
    xCy=stats::logLik(MASS::polr(yy~x,method="logistic"))
  }else{
    xCy=stats::logLik(stats::glm(yy~x,family=stats::binomial))
  }
  xCys = xCy
  improv = TRUE
  while(improv){
    improv=FALSE
    Ps=P
    for (i in 1:(nlev_y-1)){
      for (j in (i+1):nlev_y){
        tmp = Ps[i]
        Ps[i]=Ps[j]
        Ps[j]=tmp
        yy=y
        levels(yy)=as.character(Ps)
        yy=as.factor(as.numeric(as.character(yy)))
        if (nly>2){
          xCys=stats::logLik(MASS::polr(yy~x,method="logistic"))
        }else{
          xCys=stats::logLik(stats::glm(yy~x,family=stats::binomial))
        }
        if (xCys>xCy){
          P=Ps
          xCy=xCys
          improv=TRUE
        }
        tmp = Ps[i]
        Ps[i]=Ps[j]
        Ps[j]=tmp
      }
    }
  }
  return(list(M=xCy,P=P))
}

#' @title Causal Discovery for Bivariate Categorical Data
#'
#' @description Estimate the causal direction between two categorical variables.
#'
#' @param y factor, a potential effect variable
#' @param x factor, a potential cause variable
#' @param algo exhaustive search (algo="E") of category ordering or greedy search (algo="G")
#' @return A list of length 3. cd = 1 if x causes y; cd = 0 otherwise. P is the optimal ordering of the effect variable. epsilon is the difference in log-likelihood favoring x causes y.
#'
#' @export
#'
#' @examples
#' fit = COLP(CatPairs[[1]][[1]]$Diffwt,CatPairs[[1]][[1]]$Treat,algo="E")
#' fit$cd

COLP=function(y,x,algo="E"){
  if (algo=="G"){#greedy
    xCy=COLP_greedy(y,x)
    yCx=COLP_greedy(x,y)
  }else if (algo=="E"){#exhaustive
    xCy=COLP_exhaustive(y,x)
    yCx=COLP_exhaustive(x,y)
  }

  if (nlevels(x)>2){
    xCy_optim = xCy$M+stats::logLik(MASS::polr(x~1,method="logistic"))
  }else{
    xCy_optim = xCy$M+stats::logLik(stats::glm(x~1,family=stats::binomial))
  }
  if (nlevels(y)>2){
    yCx_optim = yCx$M+stats::logLik(MASS::polr(y~1,method="logistic"))
  }else{
    yCx_optim = yCx$M+stats::logLik(stats::glm(y~1,family=stats::binomial))
  }

  if (xCy_optim>yCx_optim){
    cd=1 #x -> y
    P = xCy$P
  }else if (xCy_optim<yCx_optim){
    cd=0 #y -> x
    P = yCx$P
  }else{
    cd=P=NA
  }
  return(list(epsilon=xCy_optim-yCx_optim,cd=cd,P=P))
}
