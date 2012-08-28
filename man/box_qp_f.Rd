\name{box_qp_f}
\alias{box_qp_f}
\title{
box constrained Quadratic Program (QP)
}
\description{
\code{box_qp_f} solves the minimization problem 
\deqn{\mathrm{minimize}_{u}\;\;  (b +u)' Q (b + u );\;\;  \mathrm{subject\;\; to}\;\;  \|u\|_\infty \leq \rho } 
where \eqn{Q_{m \times m}} is symmetric PSD, \eqn{u,b \in \Re^m}. The algorithm used is one-at-a-time cyclical coordinate descent.
}
\usage{
box_qp_f(Q, u, b, rho, Maxiter, tol = 10^-4)
}

\arguments{
\item{Q}{
 (Required) is a symmetric PSD matrix of dimension \eqn{m \times m}. This is a problem parameter.
}
\item{u}{
 (Required) is the optimization variable, a vector of length m. 
The value of \code{u} serves as an initialization for the coordinate-wise algorithm.

If a suitable starting point is unavailable, start with \code{u} = 0
}

\item{b}{
 (Required) is a vector of length m, this is a problem parameter. 
}

\item{rho}{
 (Required) is the degree of shrinkage. This is a non-negative scalar. 
}

\item{Maxiter}{
(Required) is an integer denoting the maximum number of iterations (full sweeps across all the m variables), to be performed by 
\code{box_qp_f}.
}
 \item{tol}{
 is the convergence tolerance. It is a real positive number (defaults to 10^-4).  
\code{box_qp_f} converges if the relative difference of the objective values is less than \code{tol}.
}
}
\details{
This box QP function is a R wrapper to a Fortran code. This is primarily meant 
to be called from the R function \code{dpglasso}.  
One needs to be very careful (as in supplying the inputs of the progra properly) while using this as a stand alone program. 
}

\value{
\item{u}{the optimal value of the argument \code{u}, upon convergence}
\item{grad_vec}{the gradient of the objective function at \code{u}}

}
\references{
This algorithm is used as a part of the algorithm DPGLASSO described in our paper:
``The Graphical Lasso: New Insights and Alternatives by Rahul Mazumder and Trevor Hastie"
available at \url{http://arxiv.org/abs/1111.5479}
}
\author{
Rahul Mazumder and Trevor Hastie
}

%\note{
%}

%\seealso{
%	    }
\examples{
set.seed(2008)

# create problem data

m<-20;
aa<-array(rnorm(m^2),dim=c(m,m));
Q<-aa%*%t(aa); 
Q<- Q + diag(rep(0.1,m));

b<-rnorm(m);

soln<-box_qp_f(Q, u=rep(0,m), b, rho=.2 , Maxiter=1000, tol = 10^-4)


}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{box constrained QP}
\keyword{coordinate descent}% __ONLY ONE__ keyword per line
