
Gem::Specification::new do |s|
  version = "1.0.0"

  files = Dir.glob("**/*") - [ 
                               Dir.glob("simple-proj-carray-*.gem"), 
                               Dir.glob("doc/**/*"),
                               Dir.glob("test/**/*"),
                             ].flatten

  s.platform    = Gem::Platform::RUBY
  s.name        = "simple-proj-carray"
  s.summary     = "An extension library for PROJ for carray"
  s.description = <<-HERE
    An extension library for PROJ for carray
  HERE
  s.version     = version
  s.licenses    = ['MIT']
  s.author      = "Hiroki Motoyoshi"
  s.email       = ""
  s.homepage    = 'https://github.com/himotoyoshi/simple-proj-carray'
  s.files       = files
  s.extensions  = [ "ext/extconf.rb" ]
  s.required_ruby_version = ">= 1.8.1"
  s.add_runtime_dependency 'simple-proj', '~> 1.0'
  s.add_runtime_dependency 'carray', '~> 1.3'
end
