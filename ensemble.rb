#!/usr/bin/env ruby

# Copyright 2016 IBM Corp.
# Author: Robert J. Prill <rjprill@us.ibm.com>
#
# Reads species abundance from species/datasets1_*.txt and writes species/dataset1_COMMUNITY.txt


class Array
  def sum
    self.inject(0) { |sum, x| sum + x }
  end
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
    end
    results << Hash[header.zip(values)]
  end
  results
end

samples = read_csv("samples.csv")
STDERR.puts samples.first

# read datasets.csv
datasets = read_csv("datasets.csv")
STDERR.puts datasets.first

levels = ["genus", "species", "subspecies"][1..1]  # LOOK: species only

levels.each do |level|
  taxnames = Hash.new
  samples.each do |s|
    sample_tag = s[:sample_tag]
    if s[:pretty_sample_tag] && !s[:pretty_sample_tag].empty?
      sample_tag = s[:pretty_sample_tag]
    end
    outfile = "output/#{level}/#{s[:group]}_#{sample_tag}_COMMUNITY.txt"
    STDERR.puts outfile
    
    community = Hash.new(0)
    datasets.each do |d|
      if d[:in_ensemble]
        i = 1
        infile = "#{s[:group]}_#{s[:sample_tag]}_#{d[:name]}.txt"
        STDERR.puts "\t" + infile
        file = "output/#{level}/#{infile}"
        if File.exists?(file)
          prediction_lines = File.readlines(file)
          previous_value = nil
          prediction_lines.each do |line|
            taxid, value, fraction, level2, name = line.chomp.split("\t")
            unless taxid.nil? || taxid.empty?  # happens when prediction file is empty
              taxnames[taxid] = name  # remember the name for later
              if previous_value && value != previous_value
                i += 1  # only increment counter if there is no tie (handle ties)
              end
              community[taxid] += 1.0 / i  # this is the rule for combining prediction lists
              previous_value = value
            end
          end
        end
      end
    end
    
    total = community.values.sum
    
    community = community.sort_by { |k, v| -v }
    
    File.open(outfile, "w") do |f|
      community.each do |taxid, value|
        pretty_value = sprintf("%0.5f", value)
        pretty_fraction = sprintf("%0.5f", value / total)
        f.puts [taxid, pretty_value, pretty_fraction, level, taxnames[taxid]].join("\t")
      end
    end
  end
end

