\name{C85}
\alias{C85}
\docType{data}
\title{
Helicase Dependent Amplification of Vimentin using the VideoScan Platform
}
\description{
A Helicase Dependent Amplification (HDA) of Vimentin (Vim) was performed. The VideoScan Platform 
(Roediger et al. (2013)) was used to monitor the amplification. The HDA was performed at 65 
degree Celsius. Three concentrations of input DNA (D1, D2, D3) were used.
}
\usage{data(C85)}
\format{
  A data frame with 301 observations on the following 5 variables.
  \describe{
    \item{\code{Cycle}}{Cycles HDA measurements.}
    \item{\code{t.D1}}{Dilution 1, elapsed time during HDA in seconds.}
    \item{\code{MFI.D1}}{Dilution 1, fluorescence.}
    \item{\code{t.D2}}{Dilution 2, elapsed time during HDA in seconds.}
    \item{\code{MFI.D2}}{Dilution 2, fluorescence.}
    \item{\code{t.D3}}{Dilution 3, elapsed time during HDA in seconds.}
    \item{\code{MFI.D3}}{Dilution 3, fluorescence.}
  }
}
\details{

To perform an isothermal amplification in VideoScan, standard conditions for the IsoAmp(R) III 
Universal tHDA Kit (Biohelix) were used. Primers and templates are described in Roediger et al. 
(2013). The reaction was composed of reaction mix A)10 micro L A. bidest, 1.25 micro L 10xbuffer, 
0.75 micro L primer(150nM final), 0.5 micro L template plasmid. Preincubation: This mixture was 
incubated for 2 min at 95 degree. Celsius and immediately placed on ice. Reaction mix B) 5 micro L 
A. bidest., 1.25 micro L 10x buffer, 2 micro L NaCl, 1.25 micro L MgSO4, 1,75 micro L dNTPs, 0.25 
micro L EvaGreen, 1 micro L enzyme mix. The mix was covered with 50 micro L mineral oil. The 
fluorescence measurement in VideoScan HCU started directly after adding buffer B at 65 degree 
Celsius. A 1x (D1), a 1:10 dilution (D2) and a 1:100 (D3) dilution were tested.
Temperature profile (after Preincubation):
- 60 seconds at 65 degree Celsius 
- 11 seconds at 55 degree Celsius && Measurement
}
\source{
Claudia Deutschmann & Stefan Roediger, BTU Cottbus - Senftenberg, Senftenberg, Germany
}
\references{
A Highly Versatile Microscope Imaging Technology Platform for the Multiplex Real-Time Detection of 
Biomolecules and Autoimmune Antibodies. S. Roediger, P. Schierack, A. Boehm, J. Nitschke, I. Berger, 
U. Froemmel, C. Schmidt, M. Ruhland, I. Schimke, D. Roggenbuck, W. Lehmann and C. Schroeder.
\emph{Advances in Biochemical Bioengineering/Biotechnology}. 133:33--74, 2013.
\url{http://www.ncbi.nlm.nih.gov/pubmed/22437246}
}
\examples{
data(C85)
# First example
plot(NA, NA, xlim = c(0,85), ylim = c(0,1), xlab = "Time [min]", 
     ylab = "Fluorescence", main = "HDA amplification")
  points(C85[, 2]/60, C85[, 3], type = "b", col = 1, pch = 20)
  points(C85[, 4]/60, C85[, 5], type = "b", col = 2, pch = 20)
  points(C85[, 6]/60, C85[, 7], type = "b", col = 3, pch = 20)
legend(40, 0.5, c("D1, 1x", "D2, 1:10", "D3, 1:100"), col = c(1:3), 
	pch = rep(20,3))

# Second example
plot(NA, NA, xlim = c(0,30), ylim = c(0,0.8), xlab = "Time [min]", 
     ylab = "Fluorescence", main = "HDA amplification")
points(C85[, 2]/60, C85[, 3], type = "b", col = 1, pch = 20)
points(C85[, 2]/60, smoother(C85[, 2]/60, C85[, 3], 
			     method = "savgol"), type = "b", col = 2, pch = 20)
points(C85[, 2]/60, smoother(C85[, 2]/60, C85[, 3], 
			     method = "smooth"), type = "b", col = 3, pch = 20)
points(C85[, 2]/60, smoother(C85[, 2]/60, C85[, 3], 
			     method = "mova"), type = "b", col = 4, pch = 20)

legend(1, 0.8, c("D1, raw", "D1, savgol", "D1, smooth", "D1, mova"), col = c(1:4), 
        pch = rep(20,4))

}
\keyword{ datasets }
\keyword{ HDA }
\keyword{ EvaGreen }
\keyword{ HCU }
\keyword{ VIM }
\keyword{ VideoScan }