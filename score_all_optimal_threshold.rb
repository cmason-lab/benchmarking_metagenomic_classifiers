#!/usr/bin/env ruby

# Author: Robert J. Prill <rjprill@us.ibm.com>
# Copyright IBM Corp 2016

# OUTPUT:
# Summary data for plotting sample HC1 SPECIES thresholding

def mkdir(dir)
  Dir.mkdir(dir) unless Dir.exist?(dir)
end
#mkdir("output")

# read samples.csv
samples = Array.new
header = nil
File.open("../samples.csv", "r").each_line do |line|
  unless header
    header = line.chomp.split(",").map { |x| x.to_sym }
    next
  end
  samples << Hash[header.zip(line.chomp.split(","))]
end
# STDERR.puts samples.first

# read datasets.csv
datasets = Array.new  
header = nil
File.open("../datasets.csv", "r").each_line do |line|
  unless header
    header = line.chomp.split(",").map { |x| x.to_sym }
    next
  end
  datasets << Hash[header.zip(line.chomp.split(","))]
end
# STDERR.puts datasets.first

levels = ["genus", "species", "subspecies"][1..1]  # LOOK: ONLY PROCESS SPECIES ONLY

levels.each do |level|
  outfile = "output/#{level}_precision_recall.txt"
  File.open(outfile, "w") do |f|

  samples.each do |s|
    if s[:is_truth_set]
      if s[:sample_tag] == "ds.soil"  # LOOK: ONLY PROCESS HC1 sample
        
        STDERR.puts s[:sample_tag]
        truth_file = "#{s[:group]}_#{s[:sample_tag]}_TRUTH.txt"
        truth = File.readlines("../output/#{level}/#{truth_file}").join
        
        if truth.split("\n").size >= 1  # more than zero entries in the truth set
          
            datasets.each do |d|
              new_name = "#{s[:group]}_#{s[:sample_tag]}_#{d[:name]}.txt"
              STDERR.puts "\t" + new_name
              file = "../output/#{level}/#{new_name}"
              if File.exists?(file)
                prediction = File.readlines(file).join
                
                # call external script to do the actual scoring
                cmd = "./score_one_optimal_threshold.rb ../output/#{level}/#{truth_file} ../output/#{level}/#{new_name}"
                lines = `#{cmd}`.chomp.split("\n")
                
                lines.each do |line|
                  f.puts [s[:sample_tag], d[:name], line].join("\t")
                end
              end
            end
          end
        end
      end
    end
  end
end
