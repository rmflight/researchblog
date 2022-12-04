library(dplyr)
library(ggplot2)
rmf_counts = readRDS(here::here("_posts", "2021-11-25-recreating-correlation-values-from-another-manuscript", "rmf_yeast_biocounts.rds"))
rmf_info = readRDS(here::here("_posts", "2021-11-25-recreating-correlation-values-from-another-manuscript", "rmf_yeast_bioreps_info.rds")) %>%
  dplyr::mutate(sample = Sample,
                rep = BiolRep,
                sample_rep = biosample)


barton_counts_info = readRDS(here::here("_posts", "2021-11-25-recreating-correlation-values-from-another-manuscript", "yeast_counts_info.rds"))
barton_counts = barton_counts_info$counts
barton_info = barton_counts_info$info

calculate_summarize_correlations = function(in_counts, in_info){
  # in_counts = barton_counts
  # in_info = barton_info

  cor_vals = list(
    raw = cor(in_counts),
    raw_no0 = {
      in_tmp = in_counts
      in_tmp[in_counts == 0] = NA
      cor(in_tmp, use = "pairwise.complete.obs")
    },
    log = cor(log1p(in_counts)),
    log_no0 = {
      in_tmp = in_counts
      in_tmp[in_counts == 0] = NA
      cor(log(in_tmp), use = "pairwise.complete.obs")
    }
  )

  med_vals = purrr::imap_dfr(cor_vals, function(in_cor, cor_id){
    # in_cor = cor_vals[[3]]
    # cor_id = "log"
    message(cor_id)

    tmp_cor = in_cor[in_info$sample_rep, in_info$sample_rep]
    med_cor = visualizationQualityControl::median_correlations(tmp_cor, in_info$sample)
    med_cor$which = cor_id
    med_cor
  })
}


rmf_medians = calculate_summarize_correlations(rmf_counts, rmf_info)
rmf_ranges = rmf_medians %>%
  dplyr::group_by(which, sample_class) %>%
  dplyr::summarise(high = max(med_cor),
                   low = min(med_cor))
barton_medians = calculate_summarize_correlations(barton_counts, barton_info)
barton_ranges = barton_medians %>%
  dplyr::group_by(which, sample_class) %>%
  dplyr::summarise(high = max(med_cor),
                   low = min(med_cor))
out_ranges = list(
  rmf_medians = rmf_medians,
  rmf_ranges = rmf_ranges,
  barton_medians = barton_medians,
  barton = barton_ranges
)

saveRDS(out_ranges, here::here("_posts", "2021-11-25-recreating-correlation-values-from-another-manuscript", "ranges.rds"))

keep_barton = barton_medians %>%
  dplyr::filter(which %in% "raw_no0") %>%
  visualizationQualityControl::determine_outliers(.) %>%
  dplyr::filter(!outlier) %>%
  dplyr::pull(sample_id)

barton_info2 = barton_info %>%
  dplyr::filter(sample_rep %in% keep_barton)
barton_counts = barton_counts[, barton_info2$sample_rep]
count_summary = visualizationQualityControl::summarize_data(t(barton_counts), barton_info2$sample)
saveRDS(count_summary, here::here("_posts", "2021-11-25-recreating-correlation-values-from-another-manuscript", "summary.rds"))

single_plot = as.data.frame(barton_counts) %>%
  ggplot(aes(x = Snf2.01, y = Snf2.02)) +
  geom_point()

saveRDS(single_plot, here::here("_posts", "2021-11-25-recreating-correlation-values-from-another-manuscript", "raw_plot.rds"))

single_log = as.data.frame(barton_counts) %>%
  ggplot(aes(x = log1p(Snf2.01), y = log1p(Snf2.02))) +
  geom_point()
saveRDS(single_log, here::here("_posts", "2021-11-25-recreating-correlation-values-from-another-manuscript", "log_plot.rds"))
