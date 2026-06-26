## This script makes a decay curve for all putative ME previously screened


putative_MEsetdt <- rbindlist(lapply(unique(mcols(putativeME_GR)$set), function(nm) {
  gr   <- putativeME_GR[mcols(putativeME_GR)$set == nm]
  hits <- findOverlaps(gr, geomMeanGR)
  data.table(
    alpha_geomean = geomMeanGR$alpha_geomean[subjectHits(hits)],
    ME            = nm
  )
}))

# add mQTL controls as negative control
mQTL_dt <- data.table(
  alpha_geomean = geomMeanGR$alpha_geomean[
    subjectHits(findOverlaps(
      makeGRfromMyCpGPos(mQTLcontrols_hg38, "mQTLcontrols"),
      geomMeanGR))],
  ME = "mQTLcontrols"
)

putative_MEsetdt <- rbind(putative_MEsetdt, mQTL_dt)
putative_MEsetdt[, ME := relevel(factor(ME), ref = "mQTLcontrols")]

# ── Plot ──────────────────────────────────────────────────────────────────────
p_decay_putative <- plot_decay_curve(
  putative_MEsetdt,
  title = "Decay curve — putative MEs vs mQTL controls"
)

p_decay_putative