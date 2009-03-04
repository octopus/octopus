#!/usr/bin/env ruby

require 'mkmf'
extension_name = 'example1'
dir_config(extension_name)

if have_header('spglib.h') and have_library("symspg")
	create_makefile(extension_name)
end
