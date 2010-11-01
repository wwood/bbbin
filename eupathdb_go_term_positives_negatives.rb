#!/usr/bin/env ruby

require 'go'
require '../forays/eupathdb/eu_path_d_b_gene_information_table'

if __FILE__ == $0
  # nucleus subsumer
  subsumer = Bio::Go.new.subsume_tester('GO:0005634')
  # nucleus parents - these annotations are ignored
  exclusion_list = %w(GO:0005575
GO:0005622
GO:0005623
GO:0005634
GO:0043226
GO:0043227
GO:0043229
GO:0043231
GO:0044424
GO:0044464)
  
  gene_info_file = EuPathDBGeneInformationTableTest
end