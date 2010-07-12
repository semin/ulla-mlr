#!/usr/bin/env ruby

require 'rubygems'
require 'bioruby'
# Hashes to store environment restraints
# key: entry ID
# value: gap included environtmental restraints annotation
left_env_seq  = {}
left_env_sa   = {}
left_env_sse  = {}
right_env_seq = {}
right_env_sa  = {}
right_env_sse = {}

# read a tem file
tem_file  = ARGV[0]
tem_obj   = Bio::FlatFile.auto(tem_file)


tem_obj.each_entry do |entry|
  if ((entry.defintion == 'sequence') || (entry.defintion == 'structure'))
    left_aa   = '*'
    right_aa  = '*'
    left_aas  = []
    right_aas = []

    seq = entry.data.data.gsub("\n", "")
    aas = seq.split("")
    aas.each_with_index do |aa, i|
      if (aa == '-')
        '-'
      else
        left
      end
    end
  elsif (entry.defintion == 'secondary structure and phi angle')
  elsif (entry.defintion == 'solvent accessibility')
  else
    next
  end
end

