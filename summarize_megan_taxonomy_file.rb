#!/usr/bin/env ruby

# by Robert Prill <rjprill@us.ibm.com>
# Copyright 2016 IBM Corp

# USAGE:  cat lca.txt | ./summarize_megan_taxonomy_file.rb

# lca.txt has lines that look like:
# ; ; [root] root; 100; [root+] cellular organisms; 100; [SuperKingdom] Bacteria; 100; [Phylum] Chloroflexi; 100; [Class] Chloroflexia; 100; [Order] Chloroflexales; 100; [Order+] Chloroflexineae; 100; [Family] Chloroflexaceae; 100; [Genus] Chloroflexus; 100; [Species] Chloroflexus aurantiacus; 50; [Species+] Chloroflexus aurantiacus J-10-fl; 50;


# ##################################################
# settings

#cutoff_fraction = 0.001  # MEGAN default -- taxon must have this fraction of reads to be reported
#cutoff_fraction = 0.0001  # more liberal cutoff

cutoff_fraction = ARGV[0].to_f

# ##################################################

mylevels = ["root", "root+", "SuperKingdom", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Species+", "Species++", "Species+++", "Species++++"]
mycodes = ["r", "i", "k", "i", "p", "c", "o", "f", "g", "s", "ss", "ss", "ss", "ss"]

taxa = Hash.new(0)

def parse(str)
  str =~ /\[(.+)\]/
  level = $1
  str =~ /\][ ]+(.+)$/
  node = $1
  [node, level]
end

STDIN.each do |line|
  fields = line.sub(/$/, " ").chomp.split("; ")[2..-1]  # throw away read identifier
  
  # build metaphlan-like name for this taxon and all parent taxa
  metaphlan_name_parts = Array.new
  
  (0..fields.size).step(2) do |i|
    node, level = parse(fields[i])
    score = fields[i+1].to_i
    
    if node
      node = node.sub(/\<.*\>/, "").sub(/ +$/, "")
    end
    
    j = mylevels.index(level)
    if j && score == 100
      metaphlan_name_parts << [mycodes[j], "__", node.gsub(/ +/, "_")].join
      metaphlan_name_substring = metaphlan_name_parts.join("|")
      taxa[metaphlan_name_substring] += 1
    end
  end
  metaphlan_name = metaphlan_name_parts.join("|")  # we don't actually use this; good for debugging
end

# sort output by abundance
taxa = taxa.sort_by { |k, v| -v }  # now an array of arrays

n_reads = taxa[0][1]

# write report to STDOUT
puts "taxon\tn_reads\tfraction_reads"  # header for output file
taxa.each do |k, v|
  if v > 0
    fraction = v.to_f / n_reads
    if fraction >= cutoff_fraction  # see comment at top of file
      puts [k, v, fraction].join("\t")
    end
  end
end

