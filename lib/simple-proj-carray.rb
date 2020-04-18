require 'simple-proj'
require 'carray'
require 'simple_proj_carray.so'

class PROJ

  def forward_carray (lon, lat, z = nil)
    if z.nil?
      argv = [lon, lat].map{|a| a.is_a?(CArray) ? a : CA_DOUBLE(a) }
      lon, lat = *argv
      ref  = argv[argv.map(&:size).to_ca.max_addr]
      x    = ref.template
      y    = ref.template
      _forward_ca_2(lon, lat, x, y)
      return x, y
    else
      argv = [lon, lat, z].map{|a| a.is_a?(CArray) ? a : CA_DOUBLE(a) }
      lat, lon, z = *argv
      ref  = argv[argv.map(&:size).to_ca.max_addr]
      x    = ref.template
      y    = ref.template
      zo   = ref.template
      _forward_ca_3(lon, lat, z, x, y, zo)
      return x, y, zo
    end
  end
  
  alias forward_lonlat_carray forward_carray
  
  def forward_latlon_carray (lat, lon, z = nil)
    return forward_carray(lon, lat, z)
  end
  
  def inverse_carray (x, y, z = nil)
    if z.nil?
      argv = [x, y].map{|a| a.is_a?(CArray) ? a : CA_DOUBLE(a) }
      x, y = *argv
      ref  = argv[argv.map(&:size).to_ca.max_addr]
      lon    = ref.template
      lat    = ref.template
      _inverse_ca_2(x, y, lon, lat)
      return lon, lat
    else
      argv = [x, y, z].map{|a| a.is_a?(CArray) ? a : CA_DOUBLE(a) }
      x, y, z = *argv
      ref  = argv[argv.map(&:size).to_ca.max_addr]
      lon  = ref.template
      lat  = ref.template
      zo   = ref.template
      _inverse_ca_3(x, y, z, lon, lat, zo)
      return lon, lat, zo
    end
  end

  alias inverse_lonlat_carray inverse_carray
  
  def inverse_lonlat_carray (lat, lon, z = nil)
    return inverse_carray(lon, lat, z)
  end

  def transform_forward_carray (x, y, z = nil)
    if z.nil?
      argv = [x, y].map{|a| a.is_a?(CArray) ? a : CA_DOUBLE(a) }
      x, y = *argv
      ref  = argv[argv.map(&:size).to_ca.max_addr]
      xo   = ref.template
      yo   = ref.template
      p x,y
      _transform_forward_ca_2(x, y, xo, yo)
      return xo, yo
    else
      argv = [x, y, z].map{|a| a.is_a?(CArray) ? a : CA_DOUBLE(a) }
      x, y, z = *argv
      ref  = argv[argv.map(&:size).to_ca.max_addr]
      xo   = ref.template
      yo   = ref.template
      zo   = ref.template
      _transform_forward_ca_3(x, y, z, xo, yo, zo)
      return xo, yo, zo
    end
  end

  def transform_inverse_carray (x, y, z = nil)
    if z.nil?
      argv = [x, y].map{|a| a.is_a?(CArray) ? a : CA_DOUBLE(a) }
      x, y = *argv
      ref  = argv[argv.map(&:size).to_ca.max_addr]
      xo  = ref.template
      yo  = ref.template
      _transform_inverse_ca_2(x, y, xo, yo)
      return xo, yo
    else
      argv = [x, y, z].map{|a| a.is_a?(CArray) ? a : CA_DOUBLE(a) }
      x, y, z = *argv
      ref  = argv[argv.map(&:size).to_ca.max_addr]
      xo   = ref.template
      yo   = ref.template
      zo   = ref.template
      _transform_inverse_ca_3(x, y, z, xo, yo, zo)
      return xo, yo, zo
    end
  end

end
