require "mkmf"
require 'carray/mkmf'

if defined?(Gem)
  if Gem::VERSION >= "1.7.0"    
    rb_proj_include = [$sitearchdir, *Gem::Specification.find_all_by_name("simple-proj").map { |s| File.join(s.gem_dir,"ext") }]
  else
    rb_proj_include = [$sitearchdir, *Gem.all_load_paths.reverse.grep(/simple-proj/)]
  end
else
  rb_proj_include = [$sitearchdir]
end

dir_config("proj", possible_includes, possible_libs)

$CFLAGS += rb_proj_include.map{|s| " -I #{s}"}.join(" ")

if have_header("rb_proj.h") 
  have_carray()
  create_makefile("simple_proj_carray")
end


