\name{dpglasso}
\alias{dpglasso}

\title{
dpglasso
}
\description{
Does block (one row/column at a time) coordinate-wise optimization on the primal of the Graphical Lasso problem:

\deqn{\min_X  -\log\det(X) + \mbox{trace}(X\Sigma) +\rho*\|X\|_1;\;\;
  \mbox{subject to }  X \succeq 0\mbox{  (A)}}{min_X log det (X) +
  trace(X Sigma) + rho |X|_1} 

}
\usage{
dpglasso(Sigma, X=NULL,invX=NULL,rho,outer.Maxiter=100,obj.seq=FALSE,outer.tol=10^-5)
}

\arguments{
  \item{Sigma}{(Required) the sample covariance matrix,  symmetric PSD with dimensions \eqn{p \times p}.}
  \item{X}{is an initialization to the precision matrix X. It must be symmetric, PD with dimensions \eqn{p \times p}.
Defaults to  \code{X= diag( 1/(rep(rho,p) + diag(Sigma)) )} }
\item{invX}{is an initialization to the covariance matrix. 
It must be symmetric with dimensions \eqn{p \times p}. 
It is not necessary for \code{invX} to be the inverse of \code{X}.
Defaults to \code{invX<- Sigma + diag(rep(rho,p))} }

\item{rho}{
 (Required)  is the amount of regularization. It is a non-negative scalar.
}

\item{outer.Maxiter}{
the maximum number of outer iterations (i.e. row/column updates) to be performed.  

\code{outer.Maxiter} defaults to 100.
}

  \item{obj.seq}{
Logical variable taking values \code{TRUE/FALSE}. If \code{obj.seq=TRUE} dpglasso 
computes the objective value after every sweep across p rows/columns. 

\code{obj.seq} defaults to \code{FALSE}

Note: Computing the objective values is \eqn{O(p^3)}, and can take quite some time depending upon the size of the problem. 
Hence, it is not recommended to compute the objective values, during the course of the algorithm. 
                }

\item{outer.tol}{
 convergence criterion. \code{outer.tol} is a non-negative scalar. 
If relative difference in the frobenius norm of the precision matrices across two successive iterations is below \code{outer.tol}, 
algorithm \code{dpglasso} converges.
}

}

\details{
\code{dpglasso} can also be used as a path algorithm ie solve problem (A) on a grid of \code{rho} values.
In that case, the estimates of the precision matrix \code{X} and covariance matrix \code{invX} obtained by solving 
(A) for a certain \code{rho}, are to be supplied as warm-starts to solve problem (A) for a smaller value of 
\code{rho}. See the example below.
}
\value{
 \item{X}{precision matrix}
  \item{invX}{covariance matrix}
  \item{time.counter.QP}{This is a three dimensional vector, representing the total time taken to solve all the QPs; uses the R function \code{proc.time()}}
}

\references{
This algorithm DPGLASSO is described in the paper:
``The Graphical Lasso: New Insights and Alternatives by Rahul Mazumder and Trevor Hastie"
available at \url{http://arxiv.org/abs/1111.5479}
}
\author{
Rahul Mazumder and Trevor Hastie
}

%\note{
%    }

%\seealso{
%}

\examples{
set.seed(2008)

# create data

n=10; p = 5; 
X<-array(rnorm(n*p),dim=c(n,p)); # data-matrix 
Sigma=cov(X); # sample covariance matrix

q<-max(abs(Sigma[row(Sigma)> col(Sigma)])); 
rho=q*0.7;
B<-dpglasso(Sigma,rho=rho,outer.Maxiter=20,outer.tol=10^-6); 
# uses the default initializations for the covariance and precision matrices


# now solve the problem for a smaller value of rho,
# using the previous solution as warm-start
rho.new=rho*.8;
B.new<-dpglasso(Sigma,X=B$X,invX=B$invX,
rho=rho.new,outer.Maxiter=20,outer.tol=10^-6);

}
\keyword{ graphical lasso}
\keyword{ sparse inverse covariance selection}% __ONLY ONE__ keyword per line
\keyword{ DPGLASSO}
\keyword{Primal Graphical Lasso}
