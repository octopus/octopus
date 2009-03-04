#!/usr/bin/env ruby

require 'example1.so'
include Example1

lattice = [[4,0,0], [0,4,0], [0,0,3]]
position = [[0,0,0],
            [0.5,0.5,0.5],
            [0.3,0.3,0],
            [0.7,0.7,0],
            [0.2,0.8,0.5],
            [0.8,0.2,0.5]]
types = [1,1,2,2,2,2]
size = 6

puts "       [axis a]^T\n"
puts "Axes = [axis b]\n"
puts "       [axis c]\n"

lattice.each do |vec|
  printf("%f %f %f\n", vec[0], vec[1], vec[2])
end

puts
puts "Atom info  (type and position)"

position.each_with_index do |vec, i|
  printf("%d: %f %f %f\n", types[i], vec[0], vec[1], vec[2])
end

puts
puts "International space group"
puts example1(size, lattice, position, types)
