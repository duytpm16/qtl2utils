
plot_int_coeffCC <- function(f_blup, m_blup, map, int_scan, add_scan, diff_scan, lodcolumn, chr, pos){

      old.par <- par(no.readonly=T)


      par(mfrow = c(3,1), mar = c(0,0,0,0))
      par(mai=c(.0,.6,.18,.1))
      plot_coefCC(f_blup, map = map, main = paste0('Chr ',chr, 'Interaction Alle Effects'), xaxt = 'n', xlab = NA, ylab = 'Female\nQTL effects')
      par(mai=c(.1,.6,0,.1))
      plot_coefCC(m_blup, map = map, main = NA, xaxt = 'n', xlab = NA, ylab = 'Male\nQTL effects')
      par(mai=c(.5,.6,0,.1))

      plot_int(add_scan = add_scan, int_scan = int_scan, diff_scan = diff_scan, map = map, lodcolumn = lodcolumn, chr = chr)

      par(old.par)
}
