#!/home/eboyle/bin/Rscript --vanilla 

library(tidyverse)
library(viridis)

args = commandArgs(trailingOnly=TRUE)

# Read in data
nuc_data = read_tsv(args[1])
count_data = read_tsv(args[2])
accession_data = read_tsv(args[3]) %>% arrange(rank)
feature_annotations = read_tsv(args[4]) %>% rename(name = row_id)
model_data = read_tsv(args[5]) # "internal_output/clip_model_coef/{libname}.{clip_sample_label}.tsv"
input_replicate_label = args[6]
clip_replicate_label = args[7]
output_stem = args[8]

root_folder = dirname(dirname(args[5])) # internal_output
threshold_scan_folder = file.path(root_folder, "threshold_scan")
dir.create(threshold_scan_folder, showWarnings = FALSE, recursive = TRUE) # "output/threshold_scan/"
tested_window_folder = file.path(root_folder, "tested_windows")
dir.create(tested_window_folder, showWarnings = FALSE, recursive = TRUE)
enriched_window_folder = file.path(root_folder, "enriched_windows")
dir.create(enriched_window_folder, showWarnings = FALSE, recursive = TRUE)
enrichment_summary_folder = file.path(root_folder, "enrichment_summaries")
dir.create(enrichment_summary_folder, showWarnings = FALSE, recursive = TRUE)
all_reads_folder = file.path(root_folder, "all_reads")
dir.create(all_reads_folder, showWarnings = FALSE, recursive = TRUE)
fig_folder = file.path(root_folder, "figures")
dir.create(file.path(fig_folder, "threshold_scan"), showWarnings = FALSE, recursive = TRUE)
enriched_window_fig_folder = file.path(fig_folder, "enriched_windows")
dir.create(enriched_window_fig_folder, showWarnings = FALSE)
dir.create(file.path(fig_folder, "all_reads"), showWarnings = FALSE)

n_bin = 10

sample_cols = names(count_data)[7:ncol(count_data)] ### added, all the columns containing counts
not_included_cols = sample_cols[(sample_cols != input_replicate_label) & (sample_cols != clip_replicate_label)] # all columns that is not input or clip

# Prepare rankings of feature and transcript types
exon_subtypes = accession_data$exon_subtype %>% unique
protein_coding_subtype = accession_data$exon_subtype[accession_data$accession == "protein_coding"] %>% head(1)
prioritized_exon_subtypes = exon_subtypes[cumsum(exon_subtypes == protein_coding_subtype) < 1]
unprioritized_exon_subtypes = exon_subtypes[cumsum(exon_subtypes == protein_coding_subtype) >= 1]
feature_plot_order = c("CDS_SOLITARY", "CDS_START","CDS_STOP","CDS","UTR5","UTR3",paste0("EXON_", prioritized_exon_subtypes),paste0("EXON_", unprioritized_exon_subtypes),"SSB_ADJ","SSB_PROX","SS3_ADJ","SS3_PROX","SS5_ADJ","SS5_PROX","PRIMIRNA","INTRON") %>% rev
transcript_plot_order = feature_annotations %>% group_by(transcript_type_top) %>% count(sort=TRUE) %>% mutate(tname = gsub("_","\n",transcript_type_top)) %>% pull(tname) 

# Compile GC content information
count_data$pct_gc = nuc_data$`8_pct_gc`
#count_gc_data = count_data[select(count_data, matches("(IP|IN)_[0-9]+$")) %>% rowSums > 0,] %>% group_by(gc_bin = cut_number(pct_gc,10)) %>% filter(.data[[clip_replicate_label]] + .data[[input_replicate_label]] > 0)
#count_gc_data = count_data[count_data[sample_cols] %>% rowSums > 0,] %>% group_by(gc_bin = cut_number(pct_gc,n_bin)) # using all sample labels
count_gc_data = count_data[count_data[c(clip_replicate_label, input_replicate_label)] %>% rowSums > 0,] %>% group_by(gc_bin = cut_number(pct_gc,n_bin)) # using only labels of use

# Calculate baseline log odds based on window properties (GC content)
processed_count_data = count_gc_data %>% rename(input = all_of(input_replicate_label), clip = all_of(clip_replicate_label)) %>% summarize(b_log_odds = mean((clip / (clip + input))) %>% (VGAM::logitlink)) %>% inner_join(select(count_gc_data %>% rename(input = all_of(input_replicate_label), clip = all_of(clip_replicate_label)), -not_included_cols),.)

# Set betabinomial overdispersion
model_overdispersion = median(model_data$rho)

# Calculate p values once per unique set of CLIP and input counts and remerge into processed data
p_data = processed_count_data %>% group_by(clip, input, gc_bin, b_log_odds) %>%
	summarize %>% mutate(d_log_odds = log((clip + VGAM::logitlink(b_log_odds, inverse = TRUE)) / (input + 1 - VGAM::logitlink(b_log_odds, inverse = TRUE))) - b_log_odds) %>%
	mutate(pvalue = pmax(1e-13, 1 - VGAM::pbetabinom(q = clip - 1, size = clip + input, prob = (VGAM::logitlink)(b_log_odds,inverse=TRUE), rho = model_overdispersion %>% (VGAM::logitlink)(inverse=TRUE) ))) %>%
	inner_join(processed_count_data,.)

p_clip = with(p_data, sum(clip) / sum(clip + input))
# Optimize coverage threshold for calling enriched windows by 20% FDR
threshold_max = p_data %>% mutate(total_counts = input + clip) %>% arrange(desc(total_counts)) %>% head(100) %>% tail(1) %>% pull(total_counts)
threshold_data = tibble(threshold = seq(2, threshold_max), n_enriched = sapply(seq(2, threshold_max), function(i) sum((p_data %>% filter(input + clip >= i) %>% pull(pvalue) %>% p.adjust(.,"fdr")) < 0.2)) )
threshold_data_path = file.path(threshold_scan_folder, paste0(output_stem, ".threshold_data.tsv"))
write_tsv(threshold_data, threshold_data_path)

optimized_threshold = threshold_data %>% arrange(n_enriched %>% desc) %>% head(1) %>% pull(threshold)

pdf(file.path(fig_folder, "threshold_scan", paste0(output_stem, '.threshold_scan.pdf')), height = 1.5, width = 1.7) # 'output/figures/threshold_scan/', 
ggplot(threshold_data %>% mutate(replicate = clip_replicate_label), aes(threshold, n_enriched)) + theme_bw(base_size = 7) +
	geom_line() + xlab("Total read threshold") + ylab("# of hits (q < 20%)") + scale_x_log10()+
	geom_vline(xintercept=optimized_threshold, color = "#f86808") + theme(aspect.ratio = 1) + facet_wrap(~replicate) 
dev.off()

q_data = p_data %>% group_by(above_threshold = input + clip >= optimized_threshold) %>% 
	mutate(qvalue = ifelse(above_threshold, p.adjust(pvalue,"fdr"), NA)) %>% ungroup

all_reads_fractions_feature_data = q_data %>% inner_join(select(feature_annotations, name, feature_type_top)) %>% 
	group_by(feature_type_top) %>% summarize(input = sum(input), clip = sum(clip)) %>% 
	mutate(feature_group = sub("_.*","", feature_type_top)) %>%
	mutate(input = input / sum(input), clip = clip / sum(clip)) %>%
	mutate(feature_type_top = factor(feature_type_top, levels = feature_plot_order %>% rev))

pdf(file.path(fig_folder, "all_reads",paste0(output_stem, '.all_reads_fractions.feature.pdf')), height = 1.8, width = 2.8) # 'output/figures/all_reads/'
all_reads_fractions_feature_data %>%
	pivot_longer(names_to = "replicate", values_to = "fraction",-c(feature_type_top, feature_group)) %>%
	mutate(replicate = factor(replicate, levels = c("input", "clip"))) %>%
	mutate(clip_replicate = clip_replicate_label) %>%
ggplot(aes(x=feature_type_top, fraction, fill = feature_group, group = replicate)) + theme_bw(base_size = 7) + 
	geom_bar(stat = "identity",position = "dodge",size = 0.2) + 
	geom_bar(data = mutate(all_reads_fractions_feature_data %>% pivot_longer(names_to = "replicate", values_to = "fraction",-c(feature_type_top, feature_group)), replicate = factor(replicate, levels = c("input", "clip")), fraction = ifelse(replicate == "clip", 0, fraction)), stat = "identity",position = "dodge",size = 0.2, fill = "#ffffff", alpha = 0.5) + 
	scale_color_manual(values = c("white","black")) +
	theme(legend.position = "none", axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank()) +
	ylab("Fraction of reads") + xlab("Type of feature") + 
	facet_wrap(~clip_replicate)
dev.off()

write_tsv(all_reads_fractions_feature_data, file.path(all_reads_folder, paste0(output_stem, ".all_reads_fractions_feature_data.tsv")))

# Aggregate all reads by feature
all_reads_odds_feature_data = p_data %>% 
	inner_join(feature_annotations) %>%
	group_by(gc_bin) %>%
	mutate(b_log_odds = log(sum(clip) / sum(input))) %>%
	group_by(b_log_odds, gc_bin, feature_type_top) %>%
	summarize(log_odds = log((sum(clip)+ VGAM::logitlink(b_log_odds[1], inverse = TRUE)) / (1 - VGAM::logitlink(b_log_odds[1], inverse = TRUE) + sum(input))), n_windows = n()) %>%
	mutate(d_log_odds = log_odds - b_log_odds) %>%
	group_by(feature_type_top) %>% 
	summarize(d_log_odds = mean(d_log_odds), b_log_odds = mean(b_log_odds), n_windows = sum(n_windows)) %>%
	mutate(replicate = clip_replicate_label, odds_ratio = exp(d_log_odds)) %>%
	mutate(feature_type_top = factor(feature_type_top, levels = feature_plot_order %>% rev))

min_odds = all_reads_odds_feature_data %>% filter(n_windows >= 2.5 / p_clip / (1 - p_clip)) %>% .$odds_ratio %>% min
max_odds = all_reads_odds_feature_data %>% filter(n_windows >= 2.5 / p_clip / (1 - p_clip)) %>% .$odds_ratio %>% max

pdf(file.path(fig_folder, "all_reads", paste0(output_stem, '.all_reads_odds.feature.pdf')), height = 2, width = 2.8) # 'output/figures/all_reads/'
all_reads_odds_feature_data %>% filter(n_windows >= 2.5 / p_clip / (1 - p_clip)) %>%
	mutate(feature_group = sub("_.*","", feature_type_top)) %>%
ggplot(aes(feature_type_top, odds_ratio, fill = feature_group)) + theme_bw(base_size = 7) + 
	theme(panel.grid.minor = element_blank()) + geom_bar(stat = "identity") + 
	xlab("Type of feature") + ylab("Odds ratio") + scale_y_log10() + scale_x_discrete(drop=FALSE) +
	theme(legend.position = "none", axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1), panel.grid.major.x = element_blank()) +
	facet_wrap(~replicate) + geom_blank(data = tibble(odds_ratio = c(1,exp(1)), feature_type_top = c("INTRON","CDS"), feature_group = "CDS")) 
dev.off()

write_tsv(all_reads_odds_feature_data,file.path(all_reads_folder, paste0(output_stem, ".all_reads_odds_feature_data.tsv"))) #"output/all_reads/", 

# Aggregate all reads by feature and transcript
all_reads_odds_transcript_data = p_data %>% 
	inner_join(feature_annotations) %>%
	group_by(gc_bin) %>%
	mutate(b_log_odds = log(sum(clip) / sum(input))) %>%
	group_by(b_log_odds, gc_bin, feature_type_top, transcript_type_top) %>%
	summarize(log_odds = log((sum(clip)+ VGAM::logitlink(b_log_odds[1], inverse = TRUE)) / (1 - VGAM::logitlink(b_log_odds[1], inverse = TRUE) + sum(input))), n_windows = n()) %>%
	mutate(d_log_odds = log_odds - b_log_odds) %>%
	group_by(feature_type_top, transcript_type_top) %>% 
	summarize(d_log_odds = mean(d_log_odds), b_log_odds = mean(b_log_odds), n_windows = sum(n_windows)) %>%
	mutate(replicate = clip_replicate_label, odds_ratio = exp(d_log_odds)) %>%
	mutate(feature_type_top = factor(feature_type_top, levels = feature_plot_order %>% rev)) %>%
	mutate(transcript_type_top = factor(gsub("_", "\n", transcript_type_top), levels = transcript_plot_order %>% rev)) 

min_odds = all_reads_odds_transcript_data %>% filter(n_windows >= 2.5 / p_clip / (1 - p_clip)) %>% .$odds_ratio %>% min
max_odds = all_reads_odds_transcript_data %>% filter(n_windows >= 2.5 / p_clip / (1 - p_clip)) %>% .$odds_ratio %>% max

pdf(file.path(fig_folder, "all_reads", paste0(output_stem, '.all_reads_odds.all_transcript_types.pdf')), height = 2.6, width = 5.2) # 'output/figures/all_reads/'
all_reads_odds_transcript_data %>% filter(n_windows >= 2.5 / p_clip / (1 - p_clip)) %>%
	mutate(replicate = clip_replicate_label) %>% 
ggplot(aes(transcript_type_top, feature_type_top, fill = odds_ratio)) + theme_bw(base_size = 6) +  #, width = 0.9 * transcript_n / length(feature_plot_order)
	geom_tile() +
	scale_fill_gradient2("Odds ratio", low = RColorBrewer::brewer.pal(5, "PiYG")[1], mid = RColorBrewer::brewer.pal(5, "PiYG")[3], high = RColorBrewer::brewer.pal(5, "PiYG")[5], trans = "log", limits = c(pmin(min_odds, 1 / max_odds), pmax(max_odds, 1 / min_odds)), breaks = c(2 * pmin(min_odds, 1 / max_odds), 1, 0.5 * pmax(max_odds, 1 / min_odds)), labels = sprintf(fmt = "%.3g", c(2 * pmin(min_odds, 1 / max_odds), 1, 0.5 * pmax(max_odds, 1 / min_odds)))) +
	theme(legend.key.size = unit(0.25, "cm"), axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1),panel.grid.minor = element_blank()) +
	facet_wrap(~replicate) +
	xlab("Type of transcript") + ylab("Type of feature") + theme(axis.text.x = element_text(lineheight = 0.8))
dev.off()

write_tsv(all_reads_odds_transcript_data %>% mutate(transcript_type_top = gsub("\n", "_", transcript_type_top)), file.path(all_reads_folder, paste0(output_stem, ".all_reads_odds_transcript_data.tsv"))) # "output/all_reads/", 

# Aggregate all reads by feature and GC content
all_reads_odds_feature_gc_data = p_data %>% 
	inner_join(feature_annotations) %>%
	group_by(gc_bin) %>%
	mutate(b_log_odds = log(sum(clip) / sum(input))) %>%
	group_by(feature_type_top, gc_bin, b_log_odds) %>%
	summarize(input = sum(input), clip = sum(clip), log_odds = log((sum(clip) + VGAM::logitlink(b_log_odds[1], inverse = TRUE)) / (1 - VGAM::logitlink(b_log_odds[1], inverse = TRUE) + sum(input))), n_windows = n()) %>%
	mutate(d_log_odds = log_odds - b_log_odds) %>%
	group_by(feature_type_top, gc_bin) %>% 
	mutate(gc_decile = gc_bin %>% as.numeric %>% as.factor, replicate = clip_replicate_label, odds_ratio = exp(d_log_odds)) %>%
	mutate(feature_type_top = factor(feature_type_top, levels = feature_plot_order %>% rev))

min_odds = all_reads_odds_feature_gc_data %>% filter(n_windows >= 2.5 / p_clip / (1 - p_clip)) %>% .$odds_ratio %>% min
max_odds = all_reads_odds_feature_gc_data %>% filter(n_windows >= 2.5 / p_clip / (1 - p_clip)) %>% .$odds_ratio %>% max

pdf(file.path(fig_folder, "all_reads", paste0(output_stem, '.all_reads_odds.feature_gc.pdf')), height = 2, width = 2.5) #'output/figures/all_reads/'
all_reads_odds_feature_gc_data %>% 
	filter(n_windows >= 2.5 / p_clip / (1 - p_clip)) %>%
ggplot(aes(gc_decile, feature_type_top, fill = odds_ratio)) + theme_bw(base_size = 7) + 
	theme(panel.grid = element_blank()) + geom_tile() + 
	scale_fill_gradient2("Odds ratio", low = RColorBrewer::brewer.pal(5, "PiYG")[1], mid = RColorBrewer::brewer.pal(5, "PiYG")[3], high = RColorBrewer::brewer.pal(5, "PiYG")[5], trans = "log", limits = c(pmin(min_odds, 1 / max_odds), pmax(max_odds, 1 / min_odds)), breaks = c(3/2 * pmin(min_odds, 1 / max_odds), 1, 2/3 * pmax(max_odds, 1 / min_odds)), labels = sprintf(fmt = "%.3g", c(3/2 * pmin(min_odds, 1 / max_odds), 1, 2/3 * pmax(max_odds, 1 / min_odds)))) +
	xlab("GC decile") + ylab("Type of feature") +
	facet_wrap(~replicate) + geom_blank(data = tibble(gc_decile = c(1,2), odds_ratio = c(1,exp(1)), feature_type_top = c("INTRON","CDS")))
dev.off()

write_tsv(all_reads_odds_feature_gc_data, file.path(all_reads_folder, paste0(output_stem, ".all_reads_odds_feature_gc_data.tsv"))) #"output/all_reads/"

# Aggregate window data by feature type
all_window_feature_data = feature_annotations %>% group_by(feature_type_top) %>% count(name = "n_windows") %>% ungroup

tested_window_feature_data = q_data %>% filter(above_threshold) %>% select(-above_threshold) %>% mutate(across(c("b_log_odds", "d_log_odds","pvalue","qvalue"), ~ sprintf(.x, fmt = "%.6g")))
write_tsv(tested_window_feature_data, file.path(tested_window_folder, paste0(output_stem, ".tested_windows.tsv.gz")))# "output/tested_windows/", 

enriched_window_data = q_data %>% filter(above_threshold, qvalue < 0.2) %>% select(-above_threshold) %>% left_join(feature_annotations) %>% arrange(qvalue)
write_tsv(enriched_window_data %>% mutate(across(c("b_log_odds", "d_log_odds","pvalue","qvalue"), ~ sprintf(.x, fmt = "%.6g"))), file.path(enriched_window_folder, paste0(output_stem, ".enriched_windows.tsv.gz"))) # "output/enriched_windows/"

if(nrow(enriched_window_data) == 0) {
	for(output_file in c(
		file.path(enriched_window_fig_folder, paste0(output_stem, '.enriched_window_coverage.pdf')),
		file.path(enriched_window_fig_folder, paste0(output_stem, '.enriched_window_rates.pdf')),
		file.path(enriched_window_fig_folder, paste0(output_stem, '.enriched_window_counts.linear.pdf')),
		file.path(enriched_window_fig_folder, paste0(output_stem, '.enriched_window_counts.log10.pdf')),
		file.path(enriched_window_fig_folder, paste0(output_stem, '.enriched_window_odds.feature.pdf')),
		file.path(enriched_window_fig_folder, paste0(output_stem, '.enriched_window_odds.all_transcript_types.pdf')),
		file.path(enriched_window_fig_folder, paste0(output_stem, '.enriched_window_odds.select_transcript_types.pdf')),
		file.path(enriched_window_fig_folder, paste0(output_stem, '.enriched_window_counts.per_gene_feature.pdf'))
	))
	{
		pdf(output_file, height = 1, width = 2)
			print(ggplot() + annotate("text", x = 1, y = 1, label = "No enriched windows") + theme_void())
		dev.off()
	}
	file.create(file.path(enrichment_summary_folder, paste0(output_stem, ".enriched_window_feature_data.tsv")))
	file.create(file.path(enrichment_summary_folder, paste0(output_stem, ".enriched_window_transcript_data.tsv")))
	file.create(file.path(enrichment_summary_folder, paste0(output_stem, ".enriched_window_gene_data.tsv")))
	quit()
}

pdf(file.path(enriched_window_folder, paste0(output_stem, '.enriched_window_coverage.pdf')), height = 2, width = 2.8)
enriched_window_data %>%
	mutate(replicate = clip_replicate_label, `IP reads per region` = clip) %>%
	mutate(feature_group = sub("_.*","", feature_type_top), feature_type_top = factor(feature_type_top, levels = feature_plot_order %>% rev)) %>%
ggplot(aes(feature_type_top, `IP reads per region`, fill = feature_group)) + theme_bw(base_size = 7) + geom_violin(scale="width", size = 0.4) + scale_y_log10() +
	theme(legend.position = "none", axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1)) +
	xlab("Type of feature") + facet_grid(~replicate)
dev.off()

enriched_window_feature_data = left_join(
		all_window_feature_data,
		tested_window_feature_data %>% inner_join(select(feature_annotations, name, feature_type_top)) %>% group_by(feature_type_top) %>% count(name = "n_tested")) %>%
	left_join(
		enriched_window_data %>% group_by(feature_type_top) %>% count(name = "n_enriched")
	) %>% replace_na(list(n_tested = 0, n_enriched = 0)) %>%
	mutate(feature_type_top = factor(feature_type_top, levels = feature_plot_order %>% rev)) %>%
	mutate(`Test rate` = n_tested / n_windows, `Positive rate` = n_enriched / n_tested, `Enrichment rate` = n_enriched / n_windows)

enriched_window_feature_data$pseudocount = enriched_window_feature_data %>% with(nrow(.) * (n_windows / (sum(n_windows))))
enriched_window_feature_data$odds_ratio = enriched_window_feature_data %>% with((n_enriched + pseudocount) / (sum((n_enriched + pseudocount)) - (n_enriched + pseudocount)) / (n_windows / (sum(n_windows) - n_windows)))

pdf(file.path(enriched_window_fig_folder, paste0(output_stem, '.enriched_window_rates.pdf')), height = 3.1, width = 2.2)
enriched_window_feature_data %>%
	mutate(feature_group = sub("_.*","", feature_type_top)) %>% 
	mutate(feature_type_top = factor(feature_type_top, levels = feature_plot_order %>% rev)) %>%
	pivot_longer(names_to="metric", values_to = "value", c("Test rate","Positive rate", "Enrichment rate")) %>%
	mutate(metric = factor(metric, levels = c("Test rate", "Positive rate", "Enrichment rate"))) %>%
	mutate(label = clip_replicate_label) %>%
ggplot(aes(feature_type_top, value, color = feature_group, group = feature_type_top)) + theme_bw(base_size = 7) + 
	geom_point(stroke = 0) + facet_grid(metric~label) +
	theme(legend.position = "none", axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1)) +
	xlab("Type of feature") + scale_y_log10(limits = c(NA,1)) + ylab("Rate")
dev.off()

pdf(file.path(enriched_window_fig_folder, paste0(output_stem, '.enriched_window_counts.linear.pdf')), height = 1.8, width = 2.2)
enriched_window_feature_data %>% mutate(replicate = clip_replicate_label) %>% 
	mutate(feature_group = sub("_.*","", feature_type_top)) %>% 
ggplot(aes(feature_type_top, n_enriched, fill = feature_group, group = feature_type_top)) + theme_bw(base_size = 7) + 
	geom_bar(stat="identity") + facet_grid(.~replicate) +
	theme(legend.position = "none", axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1)) +
	xlab("Type of feature") + ylab("# enriched")
dev.off()

pdf(file.path(enriched_window_fig_folder, paste0(output_stem, '.enriched_window_counts.log10.pdf')), height = 1.8, width = 2.2)
enriched_window_feature_data %>% mutate(replicate = clip_replicate_label) %>% 
	mutate(feature_group = sub("_.*","", feature_type_top)) %>% 
ggplot(aes(feature_type_top, n_enriched, color = feature_group, group = feature_type_top)) + theme_bw(base_size = 7) + 
	geom_point(stroke = 0,) + facet_grid(.~replicate) +
	theme(legend.position = "none", axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1)) +
	xlab("Type of feature") + ylab("# enriched") + scale_y_log10()
dev.off()

pdf(file.path(enriched_window_fig_folder, paste0(output_stem, '.enriched_window_odds.feature.pdf')), height = 1.8, width = 2.2)
enriched_window_feature_data %>% mutate(replicate = clip_replicate_label) %>% 
	filter(n_enriched > 0 | n_windows > 100) %>%
	mutate(feature_group = sub("_.*","", feature_type_top)) %>% 
ggplot(aes(feature_type_top, odds_ratio, fill = feature_group, group = feature_type_top)) + theme_bw(base_size = 7) + 
	geom_bar(stat="identity") + facet_grid(.~replicate) +
	theme(legend.position = "none", axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1)) +
	scale_y_log10() + scale_x_discrete(drop=FALSE) +
	xlab("Type of feature") + ylab("Enrichment (OR)")
dev.off()

write_tsv(enriched_window_feature_data, file.path(enrichment_summary_folder, paste0(output_stem, ".enriched_window_feature_data.tsv")))

# Aggregate window data by transcript and feature type
all_window_transcript_data = feature_annotations %>% group_by(transcript_type_top, feature_type_top) %>% count(name = "n_windows") %>% ungroup

enriched_window_transcript_data = left_join(
		all_window_transcript_data,
		tested_window_feature_data %>% inner_join(select(feature_annotations, name, transcript_type_top, feature_type_top)) %>% group_by(transcript_type_top, feature_type_top) %>% count(name = "n_tested")) %>%
		replace_na(list(n_windows = 0, n_tested = 0)) %>%
	left_join(
		enriched_window_data %>% group_by(transcript_type_top, feature_type_top) %>% count(name = "n_enriched")
	) %>% replace_na(list(n_enriched = 0)) %>%
	mutate(feature_type_top = factor(feature_type_top, levels = feature_plot_order %>% rev)) %>%
	mutate(transcript_type_top = factor(gsub("_", "\n", transcript_type_top), levels = transcript_plot_order %>% rev)) %>%
	mutate(`Test rate` = n_tested / n_windows, `Positive rate` = n_enriched / n_tested, `Enrichment rate` = n_enriched / n_windows)

enriched_window_transcript_data$pseudocount = enriched_window_transcript_data %>% with(nrow(.) * (n_windows / (sum(n_windows))))
enriched_window_transcript_data$odds_ratio = enriched_window_transcript_data %>% with((n_enriched + pseudocount) / (sum((n_enriched + pseudocount)) - (n_enriched + pseudocount)) / (n_windows / (sum(n_windows) - n_windows)))

write_tsv(enriched_window_transcript_data %>% mutate(transcript_type_top = gsub("\n", "_", transcript_type_top)), file.path(enrichment_summary_folder, paste0(output_stem, ".enriched_window_transcript_data.tsv")))

min_odds = enriched_window_transcript_data %>% filter((n_enriched > 0 & n_windows > 20) | n_windows > 100) %>% .$odds_ratio %>% min
max_odds = enriched_window_transcript_data %>% filter((n_enriched > 0 & n_windows > 20) | n_windows > 100) %>% .$odds_ratio %>% max

pdf(file.path(enriched_window_fig_folder, paste0(output_stem, '.enriched_window_odds.all_transcript_types.pdf')), height = 2.6, width = 5.2)
enriched_window_transcript_data %>% filter((n_enriched > 0 & n_windows > 20) | n_windows > 100) %>%
	mutate(replicate = clip_replicate_label) %>%
ggplot(aes(transcript_type_top, feature_type_top, fill = odds_ratio)) + theme_bw(base_size = 6) +  #, width = 0.9 * transcript_n / length(feature_plot_order)
	geom_tile() +
	scale_fill_gradient2("Odds ratio", low = RColorBrewer::brewer.pal(5, "PiYG")[1], mid = RColorBrewer::brewer.pal(5, "PiYG")[3], high = RColorBrewer::brewer.pal(5, "PiYG")[5], trans = "log", limits = c(pmin(min_odds, 1 / max_odds), pmax(max_odds, 1 / min_odds)), breaks = c(2 * pmin(min_odds, 1 / max_odds), 1, 0.5 * pmax(max_odds, 1 / min_odds)), labels = sprintf(fmt = "%.3g", c(2 * pmin(min_odds, 1 / max_odds), 1, 0.5 * pmax(max_odds, 1 / min_odds)))) +
	theme(legend.key.size = unit(0.25, "cm"), axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1),panel.grid.minor = element_blank()) +
	facet_wrap(~replicate) + 
	xlab("Type of transcript") + ylab("Type of feature") + theme(axis.text.x = element_text(lineheight = 0.8))
dev.off()

single_feature_transcript_types = enriched_window_transcript_data %>% select(feature_type_top,transcript_type_top) %>% unique %>% group_by(transcript_type_top) %>% count(sort=TRUE) %>% filter(n == 1) %>% pull(transcript_type_top) 
min_odds = enriched_window_transcript_data %>% filter((n_enriched > 0 & n_windows > 20) | n_windows > 100, ! transcript_type_top %in% single_feature_transcript_types) %>% .$odds_ratio %>% min
max_odds = enriched_window_transcript_data %>% filter((n_enriched > 0 & n_windows > 20) | n_windows > 100, ! transcript_type_top %in% single_feature_transcript_types) %>% .$odds_ratio %>% max

pdf(file.path(enriched_window_fig_folder, paste0(output_stem, '.enriched_window_odds.select_transcript_types.pdf')), height = 2.6, width = 4)
enriched_window_transcript_data %>% filter(! transcript_type_top %in% single_feature_transcript_types) %>% 
	filter((n_enriched > 0 & n_windows > 20) | n_windows > 100, !transcript_type_top %in% single_feature_transcript_types) %>%
	mutate(feature_group = sub("_.*","", feature_type_top), replicate = clip_replicate_label) %>% 
ggplot(aes(transcript_type_top, feature_type_top, fill = odds_ratio)) + theme_bw(base_size = 6) +  #, width = 0.9 * transcript_n / length(feature_plot_order)
	geom_tile() +
	scale_fill_gradient2("Odds ratio", low = RColorBrewer::brewer.pal(5, "PiYG")[1], mid = RColorBrewer::brewer.pal(5, "PiYG")[3], high = RColorBrewer::brewer.pal(5, "PiYG")[5], trans = "log", limits = c(pmin(min_odds, 1 / max_odds), pmax(max_odds, 1 / min_odds)), breaks = c(2 * pmin(min_odds, 1 / max_odds), 1, 0.5 * pmax(max_odds, 1 / min_odds)), labels = sprintf(fmt = "%.3g", c(2 * pmin(min_odds, 1 / max_odds), 1, 0.5 * pmax(max_odds, 1 / min_odds)))) +
	theme(legend.key.size = unit(0.25, "cm"), axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1),panel.grid.minor = element_blank()) +
	facet_wrap(~replicate) + 
	xlab("Type of transcript") + ylab("Type of feature") + theme(axis.text.x = element_text(lineheight = 0.8))
dev.off()

# Aggregate enriched windows by feature and gene
enriched_window_gene_data = enriched_window_data %>% group_by(gene_name,feature_type_top) %>% count(name = "n_enriched") %>% 
		group_by(feature_type_top, n_enriched) %>% mutate(n_enriched = pmin(n_enriched, 5)) %>% count(name="n_genes") %>% 
		mutate(n_enriched = ifelse(n_enriched <5, as.character(n_enriched), "5+"))

pdf(file.path(enriched_window_fig_folder,paste0(output_stem, '.enriched_window_counts.per_gene_feature.pdf')), height = 2.2, width = 2.8)
enriched_window_gene_data %>%
	mutate(feature_type_top = factor(feature_type_top, levels = feature_plot_order %>% rev)) %>%
	mutate(replicate = clip_replicate_label) %>%
ggplot(aes(feature_type_top, n_genes, fill = n_enriched)) + theme_bw(base_size = 7) + 
	geom_bar(stat = "identity", position="stack") + facet_wrap(~replicate, nrow = 1) + 
	scale_fill_viridis("# enriched\nwindows",discrete=TRUE, option="plasma") +
	ylab("# genes") + xlab("Type of feature") + 
	theme(legend.position = "top", axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank())
dev.off()

write_tsv(enriched_window_gene_data, file.path(enrichment_summary_folder, paste0(output_stem, ".enriched_window_gene_data.tsv")))

