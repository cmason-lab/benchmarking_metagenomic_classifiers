#!/usr/bin/env ruby

# Author: Robert J. Prill <rjprill@us.ibm.com>
# Copyright 2016 IBM Corp.

if ARGV.size < 2
  puts "[USAGE] #{__FILE__} [truth_set.txt] [prediction.txt]"
  exit
end

truth_set_file  = ARGV[0]
prediction_file = ARGV[1]

# load truth set
@truth = Hash.new
File.open(truth_set_file, "r").each_line do |line|
  taxid, value, fraction, level, name = line.chomp.split("\t")
  value = value.to_f
  fraction = fraction.to_f
  if value > 0
    @truth[taxid] = true
  end
end
@P = @truth.size

@tp   = 0.0
@fp   = 0.0
@k    = 0.0

def precision
  @tp / @k
end

def recall
  @tp / @P
end

def increment_counters(taxid)
  @k += 1
  if @truth[taxid]
    @tp += 1
  else
    @fp += 1
  end
end

previous_value = nil
buffer = Array.new
result = Array.new

File.open(prediction_file, "r").each_line do |line|
  taxid, value, fraction, level, name = line.chomp.split("\t")
  value = value.to_f

  if previous_value && value != previous_value
    buffer.each { |t| increment_counters(t) }
    result << [precision(), recall()].join("\t")
    buffer = []
  end

  buffer << taxid
  previous_value = value
end

unless buffer.empty?  # last one
  buffer.each { |t| increment_counters(t) }
  result << [precision(), recall()].join("\t")
end

puts result.join("\n")

