dir('Conjoint/4.ninequadrants/', '*.xlsx', full.names = T, recursive = T) %>%
  lapply(read_excel) ->
  ninequadrants

names(ninequadrants) <- rep(c('sham_vs_eHx-18', 'sham_vs_eHx-36', 'sham_vs_eHx-72', 'sham_vs_pH-18', 'sham_vs_pH-36', 'sham_vs_pH-72'), each=3)
ninequadrants %>%
  ldply(.id = 'Group') %>%
  dplyr::mutate(Cut = gsub('sham_vs_', '', Group)) %>%
  tidyr::separate(Cut, into = c('Cut', 'Time'), sep = '-') %>%
  dplyr::rename(geneLog2FC=logFC...2, metaLog2FC=logFC...5,
                geneName=geneID, metaName=metaID) ->
  ninequadrants_df

dir('Conjoint/5.diffcor/', '*.cor.xlsx', full.names = T, recursive = T) %>%
  lapply(read_excel) ->
  diffcor

names(diffcor) <- c('sham_vs_eHx-18', 'sham_vs_eHx-36', 'sham_vs_eHx-72', 'sham_vs_pH-18', 'sham_vs_pH-36', 'sham_vs_pH-72')
diffcor %>%
  ldply(.id = 'Group') %>%
  dplyr::mutate(Cut = gsub('sham_vs_', '', Group)) %>%
  tidyr::separate(Cut, into = c('Cut', 'Time'), sep = '-') ->
  diffcor_df

dir('Conjoint/2.KEGG_mapping//', '*.path.xlsx', full.names = T, recursive = T) %>%
  lapply(read_excel) ->
  mapping_path
names(mapping_path) <- c('sham_vs_eHx-18', 'sham_vs_eHx-36', 'sham_vs_eHx-72', 'sham_vs_pH-18', 'sham_vs_pH-36', 'sham_vs_pH-72')
mapping_path %>%
  ldply(.id = 'Group') %>%
  dplyr::mutate(Cut = gsub('sham_vs_', '', Group)) %>%
  tidyr::separate(Cut, into = c('Cut', 'Time'), sep = '-') ->
  mapping_path_df

mapping_path_df %>%
  tidyr::separate_rows(gene_ID, sep=';') %>%
  tidyr::separate_rows(meta_ID, sep=';') %>%
  dplyr::select(Group, '#Pathway', ko, geneName=gene_ID, metaName=meta_ID, Cut, Time) ->
  mapping_path_anno

diffcor_df %>%
  dplyr::left_join(ninequadrants_df) %>%
  dplyr::left_join(mapping_path_anno) ->
  #dplyr::filter(KEGG != '--') ->
  diffcor_Log2FC_df

diffcor_Log2FC_df %>%
  tidyr::pivot_longer(contains('Log2FC'), names_to = 'Type', values_to = 'Log2FC') %>%
  dplyr::select(Cut, Time, geneName, metaName, KEGG, Type, Log2FC, PCC, PCCP, ko, '#Pathway') %>%
  dplyr::mutate(Type=gsub('Log2FC', '', Type),
                Cor_type = ifelse(PCC > 0, '+', '-')) %>%
  dplyr::filter(ko %in% koid_both_001) ->
  diffcor_Log2FC_plot_df

diffcor_Log2FC_plot_df %>%
  mutate(Cut_Pathway = sprintf("%s (%s)", `#Pathway`, Cut)) %>%
  #dplyr::filter(abs(PCC) >= 0.85) %>%
  ggplot() +
  aes(x=PCC, y=Log2FC, shape=Type, color=Type) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, show.legend=T, linetype='dashed') +
  scale_x_continuous(breaks = c(-1.0 -0.95, -0.9, -0.85, -0.8, 0.8, 0.85, 0.9, 0.95, 1.0)) +
  facet_grid(Cut_Pathway ~ Time + Cor_type, scales = 'free_x') +
  theme(strip.text.y = element_text(angle = 0, size = rel(1.2)),
        panel.spacing.x = unit(c(1,1.5,1,1.5,1), 'lines')) ->
  diffcor_plot

ggsave('diffcor.pdf', width = 14, height = 14, device = 'pdf')
# 
# ninequadrants_df %>%
#   left_join(mapping_path_anno) %>%
#   tidyr::pivot_longer(contains('Name'), names_to = 'ID_Type', values_to = 'ID') %>%
#   tidyr::pivot_longer(contains('Log2FC'), names_to = 'Log2FC_Type', values_to = 'Log2FC') %>%
#   dplyr::mutate(ID_Type = gsub('Name', '', ID_Type),
#                 Log2FC_Type = gsub('Log2FC', '', Log2FC_Type),
#                 Cor_type = ifelse(PCC > 0, '+', '-')) %>%
#   dplyr::filter(ID_Type == Log2FC_Type) %>%
#   dplyr::filter(ko %in% koid_both_001) ->
#   ninequadrants_4parallel
#  
# ninequadrants_4parallel %>%
#   #dplyr::filter(abs(PCC) >= 0.9) %>%
#   ggplot() +
#   aes(x=Time, y=Log2FC, color=ID_Type, group = ID) +
#   geom_line() +
#   facet_grid(`#Pathway` ~ Cut + Cor_type, scales='free')
# 
# ggsave('parallel.pdf', width = 14, height = 14, device = 'pdf')
