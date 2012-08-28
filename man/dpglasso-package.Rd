\name{dpglasso-package}
\alias{dpglasso-package}
\docType{package}
\title{
Package \code{dpglasso} solves the primal formulation of the Graphical Lasso
}
\description{
This package solves the primal formulation of the Graphical Lasso problem, by solving the duals (box QPs)
of the block problems. 

The outer block-wise optimization routine is written in \code{R} with the inner QP solver written in 
Fortran with a \code{R} wrapper around it. 
}

\details{
\tabular{ll}{
Package: \tab dpglasso\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
}
There are two functions \code{dpglasso} and \code{box_qp_f}.

The function \code{dpglasso} solves the primal formulation of the Graphical Lasso. 
The function \code{box_qp_f} solves a box constrained Quadratic program via cyclical coordinate descent.
}
\author{
Rahul Mazumder and Trevor Hastie

Maintainer: Rahul Mazumder <rahul.mazumder@gmail.com>
}
\references{
\emph{The Graphical Lasso: New Insights and Alternatives} by Rahul Mazumder and Trevor Hastie, 2011,
available at \url{http://arxiv.org/abs/1111.5479}     




\emph{Sparse inverse covariance estimation with the lasso}, by J. Friedman , T. Hastie and R. Tibshirani; 2007, Biostatistics.



}
\keyword{ package }
%\seealso{
%\code{\link[http://cran.r-project.org/web/packages/glasso/index.html}{<glasso>}} 
%}

%\examples{
%}
