\encoding{utf8}
\name{ONeilLambda}

\title{
    Lambda-parameters to fit REE patterns (O'Neil 2016)
}

\description{
    Calculate the shape parameters (\lambda1 to \lambda5) that describe the shape of an REE pattern.
}

\usage{
    ONeilLambda(where=WR,verbose=F,addWR=F)
}
\arguments{
	\item{where}{A data matrix containing whole rock analyses, to be projected. By default \code{WR}}
	\item{add}{Boolean. If \code{TRUE}, the results of the calculations will be added to WR as new columns with colnames = \code{c("lambda0","lambda1","lambda2","lambda3","lambda4","Nb.REE","MSWD","La_obs/La_calc" "Ce_obs/Ce_calc","Pr_obs/Pr_calc","Nd_obs/Nd_calc","Sm_obs/Sm_calc","Eu_obs/Eu_calc","Gd_obs/Gd_calc","Tb_obs/Tb_calc","Dy_obs/Dy_calc","Ho_obs/Ho_calc","Er_obs/Er_calc","Tm_obs/Tm_calc","Yb_obs/Yb_calc","Lu_obs/Lu_calc")} }
	\item{verbose}{Boolean. If \code{TRUE}, the results will be printed to the console, otherwise they will be returned invisibly}

}

\value{
    (invisibly, unless \code{verbose = TRUE}) A matrix with one line per row in \code{where} and 21 columns (see details).
	
	It also has the side effect of (globally) affecting its results to \code{results}, where they are available to \code{\link{addResults}}, \code{\link{r2clipboard}(results)}, \code{\link{HTMLTableResults}}, etc.
}

\details{
    The "lambda parameters" of \emph{O'Neil (2016)} are geometric parameter that describe the shape of a normalized REE diagram. The original paper also offers a normalization scheme (to the C1 chondrite) that is used here.
	
	Although the details of the calculation are fairly complex, they amount to (graphically) fitting a polynomial curve to the graphical representation of a normalized REE pattern (spidergram), in the conventional log plot ordered by atomic mass (or radius). With 13 REEs (Eu is excluded from the fit, as it often shows an anomaly), up to 13 parameters could be calculated, although only the first three or four are relevant. The form of the equation tat describes a spidergram is similar to
	
	\deqn{ y(x) = \lambda_0 + \lambda_1 x + \lambda_2 x^2},
	where x is the n-th REE, and y(x) its y-position in the diagram.
	
	Thus, \lambda_0 is related to the total REE abundance, \lambda_1, a linear term, is the overall slope of the REE pattern and \lambda_2, a quadratic term, is its curvature.
	
	The following figure, reproduced from O'Neil (2016), shows the \lambda paramters (three numbers on the spidergrams) associated with REE patterns of different shapes (all having \lambda_1 = 1). The background dataset is OFB (circles) and Hawaian tholeites (crosses).
	
	\if{html}{\figure{ONeil_fig7.png}}
    \if{latex}{\figure{ONeil_fig7.png}{options: width=5.5in}}
	
	This function returns, for each sample:
	\item the 5 first \lambda terms (\code{lambda0} to \code{lambda4}) describing the shape of the REE pattern (although O'Neil indicates that \code{lambda3} and \code{lambda4} are unlikely to have a geological meaning, and in any case are representative only for high-quality analyses (better than 4 \% relative error on each value for \code{lambda3}, better than 2 \% for \code{lambda4});
	\item the number of REE present (\code{Nb.REE}), Eu excluded;
	\item the reduced \chi^2 value (\code{MSWD}), an indicator of the quality of the fit (and therefore of the presence of anomalies, other than for Eu). There is no hard and fast rule, but a MSWD > 10 probably points to at least one suspect value;
	\item the ratios between the observed and calculated REE values for each REE (\code{La_obs/La_calc}, etc.). In a perfectly smooth pattern exactly represented by a fourth-order polynomial, the observed and calculated values should be equal (i.e. \code{REE_obs/REE_calc} = 1). Note that this provides an alternative way to quantify anomalies, notably Eu: \code{Eu_obs/Eu_calc} is therefore another expression of the Eu anomaly.

    If used as a plugin, a GUI function \code{.ONeilGUI} is also supplied 
    (and linked to the menu item).
}

\references{
    O'Neil HSC, 2016. The smoothness and shapes of chondrite-normalized Rare Earth Element patterns in basalts. J.Petrol. 57(8):1463-1508.
    \href{http://dx.doi.org/10.1093/petrology/egw047}

    and in particular the spreadsheet in the electronic appendix, from which the present code was adapted (translated from VBA macro). 
}

\author{
    Jean-François Moyen \email{jfmoyen@gmail.com}
}

\examples{
    data(sazava)
	accessVar("sazava")
    ONeilLambda(verbose=T,add=T)
	binary("lambda1","lambda2")
}

\keyword{}
\concept{Calculation}
\concept{REE}
\concept{Spidegram}
\concept{Menu: Plugins}

