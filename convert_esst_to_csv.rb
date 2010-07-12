#!/usr/bin/ruby -w

require 'ostruct'
require 'facets/enumerable'

type          = :char # :num or :char something else to print categorical values using characters
matrix_file   = ARGV[0]
matrix_ext    = File.extname(matrix_file)
matrix_stem   = File.basename(matrix_file, matrix_ext)
aas           = nil
envs          = []
matrices      = []
matrix_names  = []
current_array = nil
tag           = false

IO.foreach(matrix_file) do |line|
  line.chomp!
  if line =~ /^#\s+(ACDEFGHIKLMNPQRSTV\w+)/
    # detect amino acid labels
    aas = $1.split('')
  elsif line =~ /^#\s+(.*;\w+;\w+;[T|F];[T|F])/
    # detect definitions of environmental types and their classes
    elems             = $1.split(';')
    env               = OpenStruct.new
    env.name          = elems[0]
    env.tem_values    = elems[1].split('')
    env.class_labels  = elems[2].split('')
    env.constraind   = (elems[3] == 'T' ? true : false)
    env.masked       = (elems[4] == 'T' ? true : false)
    envs << env
  elsif line =~ /^#/
    # Skip comments
    next
  elsif line =~ /^>Total/i
    # skip the last 'total' table
    break
  elsif line =~ /^>(\S+)\s+(\d+)/
    # detect a label and a number for a table
    tag = true
    current_array = []
    #current_array << aas.map { |aa| "\"#{aa}\"" }
    matrix_names  << $1
  elsif line =~ /^#{aas[-1]}\s+(.*)$/
    # detect and tag the last row of a table
    # add a row to the matrix
    tag = false
    current_array << $1.strip.split(/\s+/).map(&:to_f)
    matrices << current_array
  elsif (line =~ /^\S\s+(.*)$/) && tag
    # detect a row of a matrix
    # add a row to the matrix
    current_array << $1.strip.split(/\s+/).map(&:to_f)
  else
    raise "Something wrong!: #{line}"
  end
end

# print out column names - type 1
#colnames = []
#colnames.concat(envs.map { |e| e.name })
#colnames << "AA1" << "AA2" << "FREQ"
#puts colnames.join(", ")

## iterate matrices to print out rows
#matrices.each_with_index do |matrix, mi|
  #env_cols = matrix_names[mi].split('').map_with_index { |c, li|
    #type == :num ? envs[li].class_labels.index(c) : c
  #}
  #matrix.each_with_index do |row, ri|
    #row.each_with_index do |col, ci|
      #out = []
      #out.concat(env_cols).concat(aa_cols)
      #if type == :num
        #out << ci << ri << col
      #else
        #out << aas[ci] << aas[ri] << col
      #end
      #puts out.join(",")
    #end
  #end
#end

# print out column names - type 2
colnames = []
colnames.concat(envs.map { |e| e.name })
colnames.concat(aas.map { |aa1| aas.map { |aa2| "#{aa1}#{aa2}" } }.flatten)
puts colnames.join(", ")

# iterate matrices to print out rows
matrices.each_with_index do |matrix, mi|
  env_cols = matrix_names[mi].split('').map_with_index { |c, li|
    type == :num ? envs[li].class_labels.index(c) : c
  }
  freq_cols = []
  matrix.each_with_index do |row, ri|
    row.each_with_index do |col, ci|
      freq_cols << col
    end
  end
  out_cols  = env_cols + freq_cols
  puts out_cols.join(", ")
end
