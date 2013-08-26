#!/usr/bin/env ruby
 
$/ = ">"
ARGF.gets
while rec = ARGF.gets
  rec.chomp!
  nl = rec.index("\n")
  header = rec[0..nl-1]
  seq = rec[nl+1..-1]
  seq.gsub!(/\n/,'')
  puts [header, seq.length].join(" ")
end