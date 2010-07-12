#!/usr/bin/ruby -w

require 'matrix'
require 'ostruct'

VERSION   = "0.0.1"
esst      = ARGV[0]
csv       = ARGV[1]
csv_ext   = File.extname(csv)
csv_stem  = File.basename(csv, csv_ext)
aas       = 'ACDEFGHIKLMNPQRSTVWYJU'.split('')
pad       = 6
env_feats = []
env_index = 0

# read an original esst file and print esst header properly
IO.readlines(esst).each_with_index do |line, i|
  if line =~ /^#\s+Creator/
    puts "# Creator: ulla-mlr version #{VERSION}"
  elsif line =~ /^#\s+Creation/
    puts "# Creation Date: #{Time.now.strftime("%d/%m/%Y %H:%M")}"
  elsif line =~ /^#/
    print line
  elsif line =~ /^>/
    break
  else
    warn "Something wrong with: #{line}"
  end
end

# print newly calculated values from cvs file
IO.readlines(csv).each_with_index do |line, i|
  cols      = line.chomp.gsub('"', '').split(',').map(&:strip)
  env_id    = cols[0].to_i
  next if env_id == 0 # skip header
  env_label = cols[1]
  columns   = []
  cols[2..-1].each_slice(aas.size) { |a| columns << a }
  matrix    = Matrix.columns(columns)

  puts ">#{env_label} #{env_id-1}"
  matrix_header = ("%-3s" % '#') + aas[0..-2].inject('') { |s, a| s + ("%#{pad}s" % a) }
  puts matrix_header
  (0...matrix.row_size).each do |ri|
    matrix_row = ("%-3s" % aas[ri]) + matrix.row(ri).to_a.inject('') { |s, v|
      v.is_a?(Float) ? s + ("%#{pad}.2f" % v) : s + ("%#{pad}d" % v)
    }
    puts matrix_row
  end
end
