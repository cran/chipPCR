\name{capillaryPCR}
\alias{capillaryPCR}
\docType{data}
\title{
capillary convective PCR
}
\description{
The capillary Convective PCR (ccPCR) is a modified device of the ccPCR system proposed by Chou et 
al. 2013. 
}
\usage{data(capillaryPCR)}
\format{
  A data frame with 1844 observations on the following 10 variables.
  \describe{
    \item{\code{t.121205}}{Elapsed time during amplification}
    \item{\code{ED.121205}}{a numeric vector}
    \item{\code{t.121128}}{Elapsed time during amplification}
    \item{\code{ED.121128}}{a numeric vector}
    \item{\code{t.121130.1}}{Elapsed time during amplification}
    \item{\code{ED.121130.1}}{a numeric vector}
    \item{\code{t.121130.2}}{Elapsed time during amplification}
    \item{\code{ED.121130.2}}{a numeric vector}
    \item{\code{t.121130.3}}{Elapsed time during amplification}
    \item{\code{ED.121130.3}}{a numeric vector}
  }
}
\details{
Modified version of the capillary convective tube isothermal heater heater by Chou et al. 2011. 
As heating system a conventional block heat was used. On the top of the heating block, 
we placed for the uptake of the capillaries an aluminum block (8 mm height) in which four holes 
(3.2 mm diameter and 3.0 mm depth with round shaped bottom) were drilled. The capillaries are 
regular 100 micro L Roche LightCycler(R). These glass capillaries have a round shaped closed bottom 
(2.3 mm inner diameter and 3.2 mm outer diameter). An "ESE-Log" detector (QIAGEN Lake Constance) 
was used for the real time fluorescent measurements, which was mounted in a distance of 5-10 mm 
next to the capillary. The PCR was performed with SYBR(R) Green fluorescent intercalating dye. 
Thereof the ESE-Log has in one channel the excitation at 470 nm and the detection at 520 nm. 
The data was recorded by the FL Digital Software (QIAGEN Lake Constance) and the exported text 
based raw data.

}
\source{
Ralf Himmelreich, IMM, Mainz, Germany
}
\references{
Chou, W., Chen, P., Miao Jr, M., Kuo, L., Yeh, S. and Chen, P. (2011). Rapid DNA amplification in a 
capillary tube by natural convection with a single isothermal heater. Biotech. 50, 52-57.
}
\examples{
data(capillaryPCR)
plot(NA, NA, xlim = c(0,80), ylim = c(0,1300), xlab = "time [min]", 
     ylab = "Voltage [micro V]", main = "capillary convective PCR - Raw Data")

for (i in c(1,3,5,7)) {
  lines(capillaryPCR[, i], capillaryPCR[, i+1], type = "b", pch = 20) 
}

abline(h = 290, v = c(18, 23, 35))
legend(60,800, c("Run 1", "Run 2", "Run 3", "Control"), pch = 20, lwd = 2)
  
par(mfrow = c(2,1))

type <- c("mova", "spline", "savgol")
plot(NA, NA, xlim = c(0,80), ylim = c(-50,1100), xlab = "time [min]", 
     ylab = "Voltage [micro V]", main = "capillary convective PCR with mova, spline and savgol")
for (i in 1:3) {
  for (j in c(1,3,5,7)) {
      tmp <- data.frame(na.omit(capillaryPCR[, j]), na.omit(capillaryPCR[, j+1]))
      tmp.sm <- smoother(tmp[, 1], tmp[, 2], method = type[i])
      lines(data.frame(tmp[, 1], tmp.sm), type = "b", pch = 20, cex = 0.5, col = i)
  }
}

abline(h = 50, v = c(17.5, 21.3, 32.9))
legend(0,1000, c("mova", "spline", "savgol"), pch = 20, lwd = 2, col = c(1:3))
  
plot(NA, NA, xlim = c(10,40), ylim = c(-50,100), xlab = "time [min]", 
     ylab = "Voltage [micro V]", main = "capillary convective PCR with mova, spline and savgol")
for (i in 1:3) {
  for (j in c(1,3,5,7)) {
      tmp <- data.frame(na.omit(capillaryPCR[, j]), na.omit(capillaryPCR[, j+1]))
      tmp.sm <- smoother(tmp[, 1], tmp[, 2], method = type[i])
      lines(data.frame(tmp[, 1], tmp.sm), type = "b", pch = 20, cex = 0.5, col = i)
  }
}
abline(h = 50, v = c(17.5, 21.3, 32.9))
legend(35,80, c("mova", "spline", "savgol"), pch = 20, lwd = 2, col = c(1:3))

par(mfrow = c(1,1))

}
\keyword{ datasets }
\keyword{ capillary }
