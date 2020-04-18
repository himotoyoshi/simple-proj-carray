#include "ruby.h"
#include "rb_proj.h"
#include "carray.h"

static VALUE
rb_proj_forward_ca_2 (VALUE self, volatile VALUE vlon, volatile VALUE vlat, 
                                  volatile VALUE vx, volatile VALUE vy)
{
  Proj      *proj;
  PJ_COORD   data_in, data_out;
  
  CArray    *clon, *clat, *cx, *cy;
  double    *p1, *p2, *p3, *p4;
  ca_size_t  s1,  s2,  s3,  s4;
  int8_t    *m, *mp;
  ca_size_t  count, err_count, i;
  
  int        angular_input;

  Data_Get_Struct(self, Proj, proj);

  if ( ! proj->is_src_latlong ) {
    rb_raise(rb_eRuntimeError, "requires latlong src crs. use #transform_forward instead of #forward.");
  }

  angular_input = proj_angular_input(proj->ref, PJ_FWD);

  clon = ca_wrap_readonly(vlon, CA_DOUBLE);
  clat = ca_wrap_readonly(vlat, CA_DOUBLE);
  cx   = ca_wrap_writable(vx,   CA_DOUBLE);
  cy   = ca_wrap_writable(vy,   CA_DOUBLE);
  
  ca_attach_n(4, clon, clat, cx, cy);

  count = ca_set_iterator(4, 
                          clon, &p1, &s1,
                          clat, &p2, &s2,
                          cx,   &p3, &s3,
                          cy,   &p4, &s4);

  m = ca_allocate_mask_iterator(2, clon, clat);

  mp = m;
  err_count = 0;
  for (i=0; i<count; i++) {
    if ( *mp ) {
      err_count += 1;
      *p3 = 0.0/0.0;
      *p4 = 0.0/0.0;
    }
    else {
      if ( angular_input == 1 ) {
        data_in.lp.lam = proj_torad(*p1);
        data_in.lp.phi = proj_torad(*p2);
      }
      else {
        data_in.xy.x = (*p1);
        data_in.xy.y = (*p2);        
      }
      data_out = proj_trans(proj->ref, PJ_FWD, data_in);
      if ( data_out.xy.x == HUGE_VAL ) {
        err_count += 1;
        *p3 = 0.0/0.0;
        *p4 = 0.0/0.0;
      }
      else {
        *p3 = data_out.xy.x;
        *p4 = data_out.xy.y;          
      }
    }
    p1+=s1; p2+=s2; p3+=s3; p4+=s4; mp++;
  }

  free(m);

  ca_sync_n(2, cx, cy);
  ca_detach_n(4, clon, clat, cx, cy);

  return SIZE2NUM(err_count);
}

static VALUE
rb_proj_forward_ca_3 (VALUE self, volatile VALUE vlon, volatile VALUE vlat, volatile VALUE vzi,
                                  volatile VALUE vx, volatile VALUE vy, volatile VALUE vzo)
{
  Proj      *proj;
  PJ_COORD   data_in, data_out;
  
  CArray    *clon, *clat, *czi, *cx, *cy, *czo;
  double    *p1, *p2, *p3, *p4, *p5, *p6;
  ca_size_t  s1,  s2,  s3,  s4,  s5,  s6;
  int8_t    *m, *mp;
  ca_size_t  count, err_count, i;
  
  int        angular_input;

  Data_Get_Struct(self, Proj, proj);

  if ( ! proj->is_src_latlong ) {
    rb_raise(rb_eRuntimeError, "requires latlong src crs. use #transform_forward instead of #forward.");
  }

  angular_input = proj_angular_input(proj->ref, PJ_FWD);

  clon = ca_wrap_readonly(vlon, CA_DOUBLE);
  clat = ca_wrap_readonly(vlat, CA_DOUBLE);
  czi  = ca_wrap_readonly(vzi,  CA_DOUBLE);
  cx   = ca_wrap_writable(vx,   CA_DOUBLE);
  cy   = ca_wrap_writable(vy,   CA_DOUBLE);
  czo  = ca_wrap_writable(vzo,  CA_DOUBLE);
  
  ca_attach_n(6, clon, clat, czi, cx, cy, czo);

  count = ca_set_iterator(6, 
                          clon, &p1, &s1,
                          clat, &p2, &s2,
                          czi,  &p3, &s3,
                          cx,   &p4, &s4,
                          cy,   &p5, &s5,
                          czo,  &p6, &s6);

  m = ca_allocate_mask_iterator(3, clon, clat, czi);

  mp = m;
  err_count = 0;
  for (i=0; i<count; i++) {
    if ( *mp ) {
      err_count += 1;
      *p4 = 0.0/0.0;
      *p5 = 0.0/0.0;
      *p6 = 0.0/0.0;
    }
    else {
      if ( proj->is_src_latlong == 1 ) {
        data_in.lpz.lam = proj_torad(*p1);
        data_in.lpz.phi = proj_torad(*p2);
        data_in.lpz.z   = (*p3);
      }
      else {
        data_in.xyz.x = (*p1);
        data_in.xyz.y = (*p2);        
        data_in.xyz.z = (*p3);
      }
      data_out = proj_trans(proj->ref, PJ_FWD, data_in);
      if ( data_out.xyz.x == HUGE_VAL ) {
        err_count += 1;
        *p4 = 0.0/0.0;
        *p5 = 0.0/0.0;
        *p6 = 0.0/0.0;
      }
      else {
        *p4 = data_out.xyz.x;
        *p5 = data_out.xyz.y;
        *p6 = data_out.xyz.z;
      }
    }
    p1+=s1; p2+=s2; p3+=s3; p4+=s4; p5+=s5; p5+=s5; mp++;
  }

  free(m);

  ca_sync_n(3, cx, cy, czo);
  ca_detach_n(6, clon, clat, czi, cx, cy, czo);

  return SIZE2NUM(err_count);
}

static VALUE
rb_proj_inverse_ca_2 (VALUE self, volatile VALUE vx, volatile VALUE vy, 
                                  volatile VALUE vlon, volatile VALUE vlat)
{
  Proj      *proj;
  PJ_COORD   data_in, data_out;

  CArray    *cx, *cy, *clon, *clat;
  double    *p1, *p2, *p3, *p4;
  ca_size_t  s1,  s2,  s3,  s4;
  int8_t    *m, *mp;
  ca_size_t  count, err_count, i;

  int        angular_output;
  

  Data_Get_Struct(self, Proj, proj);

  if ( ! proj->is_src_latlong ) {
    rb_raise(rb_eRuntimeError, "requires latlong src crs. use #transform_inverse instead of #inverse.");
  }

  angular_output = proj_angular_output(proj->ref, PJ_INV);

  cx   = ca_wrap_readonly(vx,   CA_DOUBLE);
  cy   = ca_wrap_readonly(vy,   CA_DOUBLE);
  clon = ca_wrap_writable(vlon, CA_DOUBLE);
  clat = ca_wrap_writable(vlat, CA_DOUBLE);

  ca_attach_n(4, cx, cy, clon, clat);

  count = ca_set_iterator(4, 
                          cx,   &p1, &s1,
                          cy,   &p2, &s2,
                          clon, &p3, &s3,
                          clat, &p4, &s4);

  m = ca_allocate_mask_iterator(2, clon, clat);
  
  mp = m;
  err_count = 0;
  for (i=0; i<count; i++) {
    if ( *mp ) {
      err_count += 1;
      *p3 = 0.0/0.0;
      *p4 = 0.0/0.0;
    }
    else {
      data_in.xy.x = (*p1);
      data_in.xy.y = (*p2);        
      data_out = proj_trans(proj->ref, PJ_INV, data_in);
      if ( data_out.lp.lam == HUGE_VAL ) {
        err_count += 1;
        *p3 = 0.0/0.0;
        *p4 = 0.0/0.0;
      }
      else {
        if ( angular_output == 1 ) {
          *p3 = proj_todeg(data_out.lp.lam);
          *p4 = proj_todeg(data_out.lp.phi);
        }
        else {
          *p3 = data_out.xy.x;
          *p4 = data_out.xy.y;          
        }
      }
    }
    p1+=s1; p2+=s2; p3+=s3; p4+=s4;
    mp++;
  }

  free(m);

  ca_sync_n(2, clon, clat);
  ca_detach_n(4, cx, cy, clon, clat);

  return SIZE2NUM(err_count);
}


static VALUE
rb_proj_inverse_ca_3 (VALUE self, volatile VALUE vx, volatile VALUE vy, volatile VALUE vzi,
                                  volatile VALUE vlon, volatile VALUE vlat, volatile VALUE vzo)
{
  Proj      *proj;
  PJ_COORD   data_in, data_out;

  CArray    *cx, *cy, *czi, *clon, *clat, *czo;
  double    *p1, *p2, *p3, *p4, *p5, *p6;
  ca_size_t  s1,  s2,  s3,  s4,  s5,  s6;
  int8_t    *m, *mp;
  ca_size_t  count, err_count, i;

  int        angular_output;

  Data_Get_Struct(self, Proj, proj);

  if ( ! proj->is_src_latlong ) {
    rb_raise(rb_eRuntimeError, "requires latlong src crs. use #transform_inverse instead of #inverse.");
  }

  angular_output = proj_angular_output(proj->ref, PJ_INV);

  cx   = ca_wrap_readonly(vx,   CA_DOUBLE);
  cy   = ca_wrap_readonly(vy,   CA_DOUBLE);
  czi  = ca_wrap_readonly(vzi,  CA_DOUBLE);
  clon = ca_wrap_writable(vlon, CA_DOUBLE);
  clat = ca_wrap_writable(vlat, CA_DOUBLE);
  czo  = ca_wrap_writable(vzo,  CA_DOUBLE);

  ca_attach_n(6, cx, cy, czi, clon, clat, czo);

  count = ca_set_iterator(6, 
                          cx,   &p1, &s1,
                          cy,   &p2, &s2,
                          czi,  &p3, &s3,
                          clon, &p4, &s4,
                          clat, &p5, &s5,
                          czo,  &p6, &s6);

  m = ca_allocate_mask_iterator(3, clon, clat, czi);
  
  mp = m;
  err_count = 0;
  for (i=0; i<count; i++) {
    if ( *mp ) {
      err_count += 1;
      *p4 = 0.0/0.0;
      *p5 = 0.0/0.0;
      *p6 = 0.0/0.0;
    }
    else {
      data_in.xyz.x = (*p1);
      data_in.xyz.y = (*p2);
      data_in.xyz.z = (*p3);
      data_out = proj_trans(proj->ref, PJ_INV, data_in);
      if ( data_out.lpz.lam == HUGE_VAL ) {
        err_count += 1;
        *p4 = 0.0/0.0;
        *p5 = 0.0/0.0;
        *p6 = 0.0/0.0;
      }
      else {
        if ( angular_output == 1 ) {
          *p4 = proj_todeg(data_out.lpz.lam);
          *p5 = proj_todeg(data_out.lpz.phi);
          *p6 = data_out.lpz.z;
        }
        else {
          *p4 = data_out.xyz.x;
          *p5 = data_out.xyz.y;
          *p6 = data_out.xyz.z;          
        }
      }
    }
    p1+=s1; p2+=s2; p3+=s3; p4+=s4;
    mp++;
  }

  free(m);

  ca_sync_n(3, clon, clat, czo);
  ca_detach_n(6, cx, cy, czi, clon, clat, czo);

  return SIZE2NUM(err_count);
}

static VALUE
rb_proj_transform_ca_2_i (VALUE self, VALUE vx, VALUE vy,
                                      VALUE vxo, VALUE vyo, PJ_DIRECTION direction)
{
  Proj *trans;
  PJ_COORD c_in, c_out;
  CArray     *cx, *cy, *cxo, *cyo;
  double     *p1, *p2, *p4,  *p5;
  ca_size_t   s1,  s2,  s4,   s5;
  int8_t  *m, *mp;
  ca_size_t     count, err_count, i;

  Data_Get_Struct(self, Proj, trans);

  cx   = ca_wrap_readonly(vx, CA_DOUBLE);
  cy   = ca_wrap_readonly(vy, CA_DOUBLE);
  cxo  = ca_wrap_writable(vxo, CA_DOUBLE);
  cyo  = ca_wrap_writable(vyo, CA_DOUBLE);

  count = ca_get_loop_count(2, cx, cy);
  if ( count == -1 ) {
    rb_raise(rb_eRuntimeError, "invalid data_num");
  }
  
  ca_attach_n(4, cx, cy, cxo, cyo);

  ca_set_iterator(4, 
                  cx, &p1, &s1,
                  cy, &p2, &s2,
                  cxo, &p4, &s4,
                  cyo, &p5, &s5);

  m = ca_allocate_mask_iterator(2, cx, cy);

  mp = m;
  err_count = 0;
  for (i=0; i<count; i++) {
    if ( *mp ) {
      err_count += 1;
      *p4 = 0.0/0.0;
      *p5 = 0.0/0.0;
    }
    else {
      c_in.xyz.x = (*p1);
      c_in.xyz.y = (*p2);
      c_out = proj_trans(trans->ref, PJ_FWD, c_in);
      (*p4) = c_out.xyz.x;
      (*p5) = c_out.xyz.y;
      if ( c_out.xyz.x == HUGE_VAL ) {
	      err_count += 1;
      }
    }
    p1+=s1; p2+=s2;  
    p4+=s4; p5+=s5; 
    mp++;
  }

  free(m);

  ca_sync_n(2, cxo, cyo);
  ca_detach_n(4, cx, cy, cxo, cyo);

  return SIZE2NUM(err_count);
}

static VALUE
rb_proj_transform_forward_ca_2 (VALUE self, VALUE vx, VALUE vy,  
                                            VALUE vxo, VALUE vyo)
{
  return rb_proj_transform_ca_2_i(self, vx, vy, vxo, vyo, PJ_FWD);
}

static VALUE
rb_proj_transform_inverse_ca_2 (VALUE self, VALUE vx, VALUE vy, 
                                            VALUE vxo, VALUE vyo)
{
  return rb_proj_transform_ca_2_i(self, vx, vy, vxo, vyo, PJ_INV);
}

static VALUE
rb_proj_transform_ca_3_i (VALUE self, VALUE vx, VALUE vy, VALUE vz, 
                                      VALUE vxo, VALUE vyo, VALUE vzo, PJ_DIRECTION direction)
{
  Proj *trans;
  PJ_COORD c_in, c_out;
  CArray     *cx, *cy, *cz, *cxo, *cyo, *czo;
  double     *p1, *p2, *p3,  *p4,  *p5,  *p6;
  ca_size_t   s1,  s2,  s3,   s4,   s5,   s6;
  int8_t  *m, *mp;
  ca_size_t     count, err_count, i;

  Data_Get_Struct(self, Proj, trans);

  cx   = ca_wrap_readonly(vx, CA_DOUBLE);
  cy   = ca_wrap_readonly(vy, CA_DOUBLE);
  cz   = ca_wrap_readonly(vz, CA_DOUBLE);
  cxo  = ca_wrap_writable(vxo, CA_DOUBLE);
  cyo  = ca_wrap_writable(vyo, CA_DOUBLE);
  czo  = ca_wrap_writable(vzo, CA_DOUBLE);

  count = ca_get_loop_count(3, cx, cy, cz);
  if ( count == -1 ) {
    rb_raise(rb_eRuntimeError, "invalid data_num");
  }
  
  ca_attach_n(6, cx, cy, cz, cxo, cyo, czo);

  ca_set_iterator(6, 
                  cx, &p1, &s1,
                  cy, &p2, &s2,
                  cz, &p3, &s3,
                  cxo, &p4, &s4,
                  cyo, &p5, &s5,
                  czo, &p6, &s6);

  m = ca_allocate_mask_iterator(3, cx, cy, cz);

  mp = m;
  err_count = 0;
  for (i=0; i<count; i++) {
    if ( *mp ) {
      err_count += 1;
      *p4 = 0.0/0.0;
      *p5 = 0.0/0.0;
      *p6 = 0.0/0.0;
    }
    else {
      c_in.xyz.x = (*p1);
      c_in.xyz.y = (*p2);
      c_in.xyz.z = (*p3);        
      c_out = proj_trans(trans->ref, PJ_FWD, c_in);
      (*p4) = c_out.xyz.x;
      (*p5) = c_out.xyz.y;
      (*p6) = c_out.xyz.z;        
/*      if ( status < 0 ) {
	      err_count += 1;
      }
*/
    }
    p1+=s1; p2+=s2; p3+=s3; 
    p4+=s4; p5+=s5; p6+=s6;
    mp++;
  }

  free(m);

  ca_sync_n(3, cxo, cyo, czo);
  ca_detach_n(6, cx, cy, cz, cxo, cyo, czo);

  return SIZE2NUM(err_count);
}

static VALUE
rb_proj_transform_forward_ca_3 (VALUE self, VALUE vx, VALUE vy, VALUE vz, 
                                            VALUE vxo, VALUE vyo, VALUE vzo)
{
  return rb_proj_transform_ca_3_i(self, vx, vy, vz, vxo, vyo, vzo, PJ_FWD);
}

static VALUE
rb_proj_transform_inverse_ca_3 (VALUE self, VALUE vx, VALUE vy, VALUE vz, 
                                            VALUE vxo, VALUE vyo, VALUE vzo)
{
  return rb_proj_transform_ca_3_i(self, vx, vy, vz, vxo, vyo, vzo, PJ_INV);
}

void
Init_simple_proj_carray ()
{
  rb_cProj = rb_define_class("PROJ", rb_cObject);

  rb_define_private_method(rb_cProj, "_forward_ca_2", rb_proj_forward_ca_2, 4);
  rb_define_private_method(rb_cProj, "_forward_ca_3", rb_proj_forward_ca_3, 6);
  rb_define_private_method(rb_cProj, "_inverse_ca_2", rb_proj_inverse_ca_2, 4);
  rb_define_private_method(rb_cProj, "_inverse_ca_3", rb_proj_inverse_ca_3, 6);

  rb_define_private_method(rb_cProj, "_transform_forward_ca_2", rb_proj_transform_forward_ca_2, 4);
  rb_define_private_method(rb_cProj, "_transform_forward_ca_3", rb_proj_transform_forward_ca_3, 6);

  rb_define_private_method(rb_cProj, "_transform_inverse_ca_2", rb_proj_transform_inverse_ca_2, 4);
  rb_define_private_method(rb_cProj, "_transform_inverse_ca_3", rb_proj_transform_inverse_ca_3, 6);

}