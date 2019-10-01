#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]

#[macro_use]
extern crate lazy_static;

extern crate libc;
use std::sync::Mutex;

include!(concat!(env!("CARGO_MANIFEST_DIR"), "/src/netcdf_bindings.rs"));
// netcdf_const.rs contains only a definiion of netCDF variables
// It uses directs values (no 'extern static' binding)
include!(concat!(env!("CARGO_MANIFEST_DIR"), "/src/netcdf_const.rs"));

// Per the NetCDF FAQ, "THE C-BASED LIBRARIES ARE NOT THREAD-SAFE"
// So, here is our global mutex.
// Use lazy-static dependency to avoid use of static_mutex feature which 
// breaks compatibility with stable channel.
lazy_static! {
    pub static ref libnetcdf_lock: Mutex<()> = Mutex::new(());
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path;
    use std::env;
    use std::ffi;

    #[test]
    fn test_nc_open_close() {
        let mnf_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
        let test_data_path = path::Path::new(&mnf_dir).join(
            "testdata").join("simple_xy.nc");
        let f = ffi::CString::new(test_data_path.to_str().unwrap()).unwrap();
        
        let mut ncid : i32 = -999999i32;
        unsafe {
            let _g = libnetcdf_lock.lock().unwrap();
            let err = nc_open(f.as_ptr(), NC_NOWRITE, &mut ncid);
            assert_eq!(err, NC_NOERR);
            let err = nc_close(ncid);
            assert_eq!(err, NC_NOERR);
        }
    }

    #[test]
    fn test_inq_varid() {
        let mnf_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
        let test_data_path = path::Path::new(&mnf_dir).join(
            "testdata").join("simple_xy.nc");
        let f = ffi::CString::new(test_data_path.to_str().unwrap()).unwrap();
        let varname = ffi::CString::new("data").unwrap();
        
        let mut ncid : i32 = -999999i32;
        let mut varid : i32 = -999999i32;
        let mut nvars : i32 = -999999i32;
        unsafe {
            let _g = libnetcdf_lock.lock().unwrap();
            let err = nc_open(f.as_ptr(), NC_NOWRITE, &mut ncid);
            assert_eq!(err, NC_NOERR);
            let err = nc_inq_nvars(ncid, &mut nvars);
            assert_eq!(err, NC_NOERR);
            assert_eq!(nvars, 1);
            let err = nc_inq_varid(ncid, varname.as_ptr(), &mut varid);
            assert_eq!(err, NC_NOERR);
            let err = nc_close(ncid);
            assert_eq!(err, NC_NOERR);
        }
    }

    #[test]
    fn test_get_var() {
        let mnf_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
        let test_data_path = path::Path::new(&mnf_dir).join(
            "testdata").join("simple_xy.nc");
        let f = ffi::CString::new(test_data_path.to_str().unwrap()).unwrap();
        let varname = ffi::CString::new("data").unwrap();
        
        let mut ncid : i32 = -999999i32;
        let mut varid : i32 = -999999i32;
        let mut buf : Vec<i32> = Vec::with_capacity(6*12);
        unsafe {
            buf.set_len(6*12);

            let _g = libnetcdf_lock.lock().unwrap();
            let err = nc_open(f.as_ptr(), NC_NOWRITE, &mut ncid);
            assert_eq!(err, NC_NOERR);

            let err = nc_inq_varid(ncid, varname.as_ptr(), &mut varid);
            assert_eq!(err, NC_NOERR);

            let err = nc_get_var_int(ncid, varid, buf.as_mut_ptr());
            assert_eq!(err, NC_NOERR);

            let err = nc_close(ncid);
            assert_eq!(err, NC_NOERR);
        }

        for x in 0..(6*12) {
            assert_eq!(buf[x], x as i32);
        }
    }

}
