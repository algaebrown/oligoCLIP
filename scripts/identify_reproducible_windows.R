library(tidyverse)



args = commandArgs(trailingOnly=TRUE)
data_directory = args[1] # internal_output/enriched_windows/
prefix = args[2] # PUM2
if(length(args) > 2) {
	blacklist = read_tsv(args[3], col_names = c("chr","start","end","name","score","strand"), col_types = "cddcdc")
} else {
	blacklist = tibble(chr=character(),start=numeric(),end=numeric(),name=character(),score=numeric(),strand=character())
}
enriched_window_files = list.files(path = data_directory, pattern = paste0(".", prefix, ".*enriched_windows.tsv.gz"), full.names = TRUE) # avoid partial gene name matching such as TIAL1 and TIAL

qvalue_thres = 0.2
# create folders
root_folder = dirname(data_directory)
# dir.create("output/figures/reproducible_enriched_windows/", showWarnings = FALSE, recursive = TRUE)
fig_folder = file.path(root_folder, "figures","threshold_scan")
dir.create(fig_folder, showWarnings = FALSE, recursive = TRUE) 
# dir.create("output/reproducible_enriched_windows/", showWarnings = FALSE, recursive = TRUE)
repwindow_folder = file.path(root_folder, "reproducible_enriched_windows")
dir.create(repwindow_folder, showWarnings = FALSE, recursive = TRUE) 

# specify col_type to handle case of no enriched windows
enriched_window_data = enriched_window_files %>%
    setNames(sub("\\.enriched_windows\\.tsv.gz", "", basename(.))) %>% #make a df with fname, name where name = 'multiplex_HEK293*.PUM2'
    map(function(x) read_tsv(x, col_types="cddcdcdddcddddcddccccccccc") %>% mutate(name = as.character(name)) %>% anti_join(blacklist %>% select(-name))) %>% 
    Filter(function(x) nrow(x) > 0, .) %>%
    bind_rows(.id = "clip_replicate_label") 

if(nrow(enriched_window_data %>% group_by(name) %>% filter(n() > 1)) == 0) {
	pdf(paste0(fig_folder, prefix, '.reproducible_enriched_window_counts.linear.pdf'), height = 1, width = 2)
		print(ggplot() + annotate("text", x = 1, y = 1, label = "No data") + theme_void())
	dev.off()
	pdf(paste0(fig_folder, prefix, '.reproducible_enriched_window_counts.log10.pdf'), height = 1, width = 2)
		print(ggplot() + annotate("text", x = 1, y = 1, label = "No data") + theme_void())
	dev.off()
	reproducible_enriched_window_data = enriched_window_data %>% group_by(across(-c(clip_replicate_label,b_log_odds,input,clip,d_log_odds,pvalue,qvalue))) %>% summarize %>% head(0)
	write_tsv(reproducible_enriched_window_data, paste0("output/reproducible_enriched_windows/", prefix, ".reproducible_enriched_windows.tsv.gz"))
	quit()
}	
reproducible_enriched_window_data = enriched_window_data %>% group_by(across(-c(clip_replicate_label,b_log_odds,input,clip,d_log_odds,pvalue,qvalue))) %>%
	summarize(
		input_sum = sum(input),
		clip_sum = sum(clip),
		enrichment_n = sum(qvalue < qvalue_thres),
		d_log_odds_min = min(d_log_odds),
		d_log_odds_mean = mean(d_log_odds),
		d_log_odds_max = max(d_log_odds),
		p_max = max(pvalue),
		p_min = min(pvalue),
		q_max = max(qvalue),
		q_min = min(qvalue)
	) %>% filter(enrichment_n > 1) %>% arrange(d_log_odds_mean %>% desc)

pdf(paste0(fig_folder, prefix, '.reproducible_enriched_window_counts.linear.pdf'), height = 1.8, width = 2.2)
reproducible_enriched_window_data %>% 
	mutate(feature_group = sub("_.*","", feature_type_top)) %>% 
ggplot(aes(feature_type_top, fill = feature_group, group = feature_type_top)) + theme_bw(base_size = 7) + 
	geom_bar() + 
	theme(legend.position = "none", axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1)) +
	xlab("Type of feature") + ylab("# enriched")
dev.off()

pdf(paste0(fig_folder, prefix, '.reproducible_enriched_window_counts.log10.pdf'), height = 1.8, width = 2.2)
reproducible_enriched_window_data %>% 
	mutate(feature_group = sub("_.*","", feature_type_top)) %>% group_by(feature_group, feature_type_top) %>% count %>%
ggplot(aes(feature_type_top, n, color = feature_group, group = feature_type_top)) + theme_bw(base_size = 7) + 
	geom_point(stroke = 0) + 
	theme(legend.position = "none", axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1)) +
	xlab("Type of feature") + ylab("# enriched") + scale_y_log10()
dev.off()

write_tsv(reproducible_enriched_window_data, paste0(repwindow_folder, prefix, ".reproducible_enriched_windows.tsv.gz"))