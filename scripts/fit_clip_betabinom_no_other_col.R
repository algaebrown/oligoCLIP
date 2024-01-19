
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

nuc_data = read_tsv(args[1])
count_data = read_tsv(args[2])
experiment = args[3]
given_clip_replicate = args[4]
outf = args[5]
print(args[2])
print(nrow(count_data))

# redirect all outputs

coef_dir = dirname(outf)
distribution_dir = file.path(dirname(coef_dir), "figures", "clip_distributions")
dir.create(coef_dir, showWarnings = FALSE, recursive = TRUE) # "output/clip_model_coef/" or "internal_output/clip_model_coef/
dir.create(distribution_dir, showWarnings = FALSE, recursive = TRUE) # "output/figures/clip_distributions/" or "internal_output/figures/clip_distributions/


n_bin = 10

sample_cols = names(count_data)[7:ncol(count_data)]
other_clip_replicates = setdiff(sample_cols, given_clip_replicate)
write('sample columns', stdout())
write(sample_cols, stdout())
write('other clip replicates', stdout())
write(other_clip_replicates, stdout())

count_data$pct_gc = nuc_data$`8_pct_gc`

print(head(count_data))

write('nrows filter', stdout())
write(nrow(count_data[count_data[sample_cols] %>% rowSums > 0,]))
count_gc_data = count_data[count_data[sample_cols] %>% rowSums > 0,] %>% group_by(gc_bin = cut_number(pct_gc,n_bin))

clip_betabinom_fit_data = lapply(other_clip_replicates, function(other_clip_replicate)
	{
		clip_gc_data = count_gc_data[as.logical(count_gc_data[,given_clip_replicate] + count_gc_data[,other_clip_replicate] > 0),]
		processed_clip_data = clip_gc_data %>% summarize(given_clip_fraction = mean(.data[[given_clip_replicate]] / (.data[[other_clip_replicate]] + .data[[given_clip_replicate]])) %>% (VGAM::logitlink)) %>% inner_join(clip_gc_data)
		betabinom_fit = VGAM::vglm(cbind(processed_clip_data[[given_clip_replicate]], processed_clip_data[[other_clip_replicate]]) ~ processed_clip_data$given_clip_fraction, VGAM::betabinomial, trace = TRUE)
		betabinom_coefs = betabinom_fit %>% (VGAM::Coef)
		global_given_clip_fraction = sum(clip_gc_data[[given_clip_replicate]]) / sum(clip_gc_data[[other_clip_replicate]] + clip_gc_data[[given_clip_replicate]])
		distribution_data = processed_clip_data %>% group_by(clip_total = .data[[given_clip_replicate]] + .data[[other_clip_replicate]], .data[[given_clip_replicate]], given_clip_fraction, gc_bin) %>%
			count(name = "count") %>% group_by(clip_total, gc_bin) %>% mutate(pdf = count / sum(count)) %>% rowwise %>% 
			mutate(binomial = dbinom(x = .data[[given_clip_replicate]], size = clip_total, prob = global_given_clip_fraction),
				betabinomial = VGAM::dbetabinom(x = .data[[given_clip_replicate]], size = clip_total, prob = (betabinom_coefs[1] + given_clip_fraction * betabinom_coefs[3]) %>% (VGAM::logitlink)(inverse=TRUE), rho = betabinom_coefs[2] %>% (VGAM::logitlink)(inverse=TRUE) )) 
		pdf(paste0(distribution_dir, experiment, ".", given_clip_replicate, ".", other_clip_replicate, '.clip_distribution.pdf'),height = 3.5, width = 6)
		print(
			ggplot(distribution_data %>% filter(clip_total <=6) %>% 
				pivot_longer(values_to="expected",names_to="distribution",-clip_total:-pdf) %>% 
				mutate(gc_bin = as.numeric(gc_bin) * n_bin), 
				aes(.data[[given_clip_replicate]] / (clip_total), pdf)) +
			geom_area(data=distribution_data %>% filter(clip_total <=6) %>% mutate(gc_bin = as.numeric(gc_bin) * n_bin)) +
			theme_bw(base_size = 7)+
			geom_line(aes(y = expected,color = distribution),linetype="21") + 
			facet_grid(clip_total~gc_bin) + xlab(paste0("Fraction of clip reads: ", given_clip_replicate)) + ylab("PDF") + 
			scale_x_continuous(breaks = c(0,1)) + ggtitle("Binned by GC percentile")
		)
		dev.off()

		betabinom_coefs %>% as_tibble(rownames="coef") %>% transmute(coef = c("mu","rho","gc_bias_term"),value) %>% pivot_wider(names_from=coef,values_from=value) %>% mutate(global_given_replicate_fraction = global_given_clip_fraction, given_replicate = given_clip_replicate, other_replicate = other_clip_replicate, experiment = experiment)
	}
) %>% bind_rows

write_tsv(clip_betabinom_fit_data, outf)
