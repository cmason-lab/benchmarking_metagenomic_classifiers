#!/usr/bin/env ruby

# Author: Robert J. Prill <rjprill@us.ibm.com>
# Copyright 2016 IBM Corp.

def mkdir(dir)
  Dir.mkdir(dir) unless Dir.exist?(dir)
end

def read_csv(filename)
  results = Array.new
  header = nil
  File.open(filename, "r").each_line do |line|
    unless header
      header = line.chomp.split(",").map { |x| x.to_sym }
      next
    end
    values = line.chomp.split(",")
    values.each_with_index do |v, i|
      if v.downcase == "true" || v.downcase == "false"
        values[i] = v == "true"  # convert "true" string and "false" string to boolean
      end
      if v.empty?
        values[i] = nil  # convert empty string to nil
      end
    end
    results << Hash[header.zip(values)]
  end
  results
end

mkdir("output")

samples = read_csv("../samples.csv")
# STDERR.puts samples.first

datasets = read_csv("../datasets.csv")
# STDERR.puts datasets.first

# append to datasets an entry for the COMMUNITY prediction
datasets << {:id=>"99", :name=>"COMMUNITY"}

levels = ["genus", "species", "subspecies"]

levels.each do |level|
  outfile = "output/#{level}_precision_recall.txt"
  File.open(outfile, "w") do |f|
    f.puts ["sample_tag", "algorithm", "", "precision", "recall"].join("\t") 
    samples.each do |s|
      score_me = false
      if level == "subspecies"  # some samples are to be scored at subspecies level
        if s[:is_strain_truth_set]
          score_me = true
        end
      else  # some samples are the be scored at the genus and species levels
        if s[:is_truth_set]
          score_me = true
        end
      end
      
      if score_me
        STDERR.puts s[:sample_tag]
        truth_file = "#{s[:group]}_#{s[:sample_tag]}_TRUTH.txt"
        if s[:pretty_sample_tag]
          truth_file = "#{s[:group]}_#{s[:pretty_sample_tag]}_TRUTH.txt"
        end
        
        truth = File.readlines("../output/#{level}/#{truth_file}").join
        if truth.split("\n").size >= 1  # more than zero entries in the truth set
          
          datasets.each do |d|
            prediction_file = "#{s[:group]}_#{s[:sample_tag]}_#{d[:name]}.txt"
            STDERR.puts "\t" + prediction_file
            file = "../output/#{level}/#{prediction_file}"
            if File.exists?(file)
              # call external script to do the actual scoring
              cmd = "./score_one_precision_recall.rb ../output/#{level}/#{truth_file} ../output/#{level}/#{prediction_file}"
              lines = `#{cmd}`.chomp.split("\n")
              lines.each do |line|
                f.puts [s[:sample_tag], d[:name], "\t"].join("\t") + line
              end
            end
          end
        end
      end
    end
  end
end

