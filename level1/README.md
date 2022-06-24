# Processing of STEREO and (potentially) LASCO data

Up until now (June, 2022) I have been working with L0 data from LASCO and STEREO as returned by sunpy (presumably via VSO) or as retrieved from the online data archives.

This was fine to start off with, but it is not very representative of the L3 processing we will be doing for CCOR, which will start from L2 data.  The L0 LASCO and STEREO data is in unsigned integer units (DN), they are not normalized by exposure time, and they have many defects and corrupted frames.  So, movies are choppy and some processing methods like noisegate don't really work.  Furthermore, the STEREO beacon images are only 512 x 512.  If we're going to test the full CCOR L3 processing, it would be good to start with 2k x 2k images.

But there is not much L1 or L2 data available online for LASCO and STEREO.  "Like SOHO, the philosophy for SECCHI has always been for the users to apply SECCHI_PREP to their data."  

The NRL team does provide some data which are processed beyond L0, but only for the COR2 telescope, and only for three-image polarization sequences, not the total brightness images.  These images can be found in

https://stereo-ssc.nascom.nasa.gov/data/ins_data/secchi/P0/

So, that is an option.  But, another option is to run the L0-L1 processing directly on tB images using `secchi_prep.pro`, which is what the code in this directory does.   Let's see if they make better movies.  I'll start with STEREO because I already know how to run the `secchi_prep.pro' SolarSoft routine to do the L0 -> L1 processing.  If that works, I may add LASCO processing to this directory in the future.

So, I started with a command like this in sunpy to grab data for a month:

```python
qr = Fido.search(timerange, a.Source(source), instrument, detector)
```

Adding a `a.Level.one` or `a.Level.two` did not help - this returns over 4000 files, all of which are L0 as far as I can tell.  I downloaded them and removed the polarized (`*_n4*`) and beacon (`*0800*`) images.  The next steps are as follows.

1. Use `process_dir.py` to process the L0 images to L1 using the SolarSoft IDL routine `secchi_prep.pro`.  In the process, make sure that they are organized properly to be further processed by CIMP, meaning that they are all of the same resolution and files of a given day are each in a separate directory.

2. Run `ComputeBackground.py` to compute the bacground

3. Feed this into L3 CIMP processing and see if it makes better movies.



