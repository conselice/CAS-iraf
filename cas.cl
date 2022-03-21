# CAS Program version -1 by C.Conselice to compute radii, asymmetry and concentration indexes, fluxes.  PLEASE DO NOT DISTRIBUTE.
#
#  Contact conselice@gmail.com with any questions/problems
#
#  must use with IRAF and must type digiphot and then phot before running
#
#The input file has the format (for this CAS program):
#
# name, x-center, y-center, x1, x2, y1, y2 
#
# where x1,x2,y1,y2 defines the x,y region of the background box 
#
# x-center,y-center are the centers of the galaxies under study, background
# is the background level of the image.
#
# Note the following conditions on input values:
#       x1 - del > 0
#       y1 - del > 0
#       x2 + del < xdim
#       y2 + del < ydim
# If the input radius is plus coordinates exceed image size then radius
# is rest to maximum size that can fit on image.
#
# Output is:
#
# image, C, E(C), A, E(A), S, E(S), R_petr, R_half, flux, flag
# 
# image is the name of the image, C is the concentration, E(C) is the error
# A is the asymmetry value, E(A) is its error
# S is the clumpiness value, E(S) is its error
# R_petr and R_half are two radii (Petrosian Radii, Half-light radii)
# flag reveals whether the edge was hit during the CAS computation (> 0)
#
procedure cas ( imgfile, outfile )

string imgfile = ""	{ prompt = "Input file: image x y x1, x2, y1, y2" }
string outfile = ""	{ prompt = "Output file: " }

real search 	{ prompt = "Pixel step-size for minimization search" }
int dc 		{ prompt = "Maximum expected center wander for minimization" }

bool vbo1 = no	{ prompt = "Level 1 diag: obj. and calc. step" }
bool vbo2 = no	{ prompt = "Level 2 diag: minimization iterations" }

struct *mylist1		{ prompt = "wiz" }
struct *mylist2		{ prompt = "bang" }
struct *mylist3         { prompt = "zoom" }

begin

	int ix1, ix2, iy1, iy2
	int ii, jj, flag, red, nrad, rad
	string image, timage, imaget, sub, st, junk1, junk2, junk3
	string imsqr
	string sqr, ifile, ofile
        string interactive
	string rot, section
	string z
	real xt, yt, inten, dis, fdis, tdis, tabflux, c20
        real bigrad, ll, se, xo, yo, bsize
	real  t, bg, csb, isb, osb, size, otflux, ne
	real flux, kk, id, otsb, econc
	real npix, tsb, nipix, tflux, halfr, si, xi, yi
	real irad, orad, conc, name, oname, newrad
	string  cutout, ssi, tmpdat
	real symm, top, bottom, tsymm,gsymm, centr, xphot, yphot, radd
	real x, y, bgx, bgy, xp, yp, rad0, e, be, jgv
	real xcen, ycen, x1, x2, y1, y2, xorg, yorg, stddev
	real btakeout, takeout, stdd, istd, error, b
	bool cflag, iflag, oflag, halfcon, dflag
	int ixx1,ixx2,iyy1,iyy2
	int xdim,ydim,ccount,bcount
	string psect

	int dd
	string bsect
	real sm, s
        real eta8, eta5, eta2, grw2, grw5, grw8, eta2t2, grw8t2, ttakeout

	ifile = imgfile
	ofile = outfile
	mylist1 = ifile

	cutout = mktemp ( "CuT" ) // ".fits"
	ssi = mktemp ( "SsI" ) // ".fits"
	rot = mktemp ( "RoT" ) // ".fits"
	sub = mktemp ( "SuB" ) // ".fits"
	sqr = mktemp ( "SqR" ) // ".fits"
	imsqr = mktemp ( "ISqR" ) // ".fits"
	tmpdat = mktemp ( "TdA" )

	imdelete ( "cutout.fits", verify-, >& "dev$null" )
	imdelete ( cutout, verify-, >& "dev$null" )
	imdelete ( ssi, verify-, >& "dev$null" )
	imdelete ( rot, verify-, >& "dev$null" )
	imdelete ( sub, verify-, >& "dev$null" )
	imdelete ( sqr, verify-, >& "dev$null" )
	imdelete ( imsqr, verify-, >& "dev$null" )
	imdelete ( "c1.fits", verify-, >& "dev$null" )
	delete ( tmpdat, verify-, >& "dev$null" )

	  imdelete ( "CuT*fits*", verify-, >& "dev$null" )
	  imdelete ( "RoT*fits*", verify-, >& "dev$null" )
	  imdelete ( "SsI*fits*", verify-, >& "dev$null" )
	  imdelete ( "SuB*fits*", verify-, >& "dev$null" )
	  imdelete ( "SqR*fits*", verify-, >& "dev$null" )
        image = "test"

        noao
        digiphot
        apphot
        artdata
        tables
        ttools
#Unlearn packages
        unlearn imstat
        unlearn phot
        unlearn tprint

#Hard-wire some parameters:
          phot.interactive = no
          phot.verify = no
          datapars.scale = 1

        cache phot
        cache imstat
        cache tprint
        cache tae

	flag = 0
        search = 0.5
        dc = 5
        halfr = 0
        tsb = 0

        imdel ("test.fits", yes, verify=no, default_acti=yes)

        print ( '#Name           x      y       C     E(C)      A      E(A)     S     E(S)    Petr   HalfR   Flux   f', >> ofile )

#	BEGIN LOOKING AT THE OBJECTS

        mylist1 = ifile


	while (fscan(mylist1, image, x,y,x1,x2,y1,y2) !=EOF)
        {


	imdelete ( "im.fits", verify-, >& "dev$null" )
	imdelete ( "test.fits", verify-, >& "dev$null" )

	xo = x
	yo = y

	  imget ( image, "i_naxis1" )
	  xdim = int(imget.value)
	  imget ( image, "i_naxis2" )
	  ydim = int(imget.value)

	  flprc; flprc

# get image dims

          size = 500

#print (size)

############################################################################
# Background calculation
############################################################################
# Get the RMS of the Background, as well as the Flux from Background.
# This section "minimizes" the background asymmetry in the same way
# as for the source, and uses this to correct the source
# asymmetry measurement.
#
# Note that there is an assumption here: the center will not wander
# by more than dc pixels from the initial center.
# Also note that the cutout does not chance even though the center does.
# This is needed for convergence.

# define extraction box and extract cutout

	  ix1 = int ( x1)
          ix2 = int ( x2)
          iy1 = int ( y1)
          iy2 = int ( y2)

          if (ix1 - dc < 0) {ix1 = 1 + dc}
          if (ix2 - dc < 0) {ix2 = 3 + dc}
          if (iy1 - dc < 0) {iy1 = 1 + dc}
          if (iy2 - dc < 0) {iy2 = 3 + dc}

	  psect = "["//ix1-dc//":"//ix2+dc//","//iy1-dc//":"//iy2+dc//"]"
	  section = "["//dc//":"//ix2-ix1+dc//","//dc//":"//iy2-iy1+dc//"]"
	  bsect = "["//ix1//":"//ix2//","//iy1//":"//iy2//"]"

          imcopy ( image//psect, cutout, verbose- )

# get BKGND mean (btakeout) from "cutout" image section

          xphot = dc + (ix2 - ix1)/2 
          yphot = dc + (iy2 - iy1)/2 

           delete ( "cntr", verify-, >& "dev$null" )
           delete ( "cord.out", verify-, >& "dev$null" )
           delete ( "f", verify-, >& "dev$null" )
           delete ( "ff.tab", verify-, >& "dev$null" )
           delete ( "fff", verify-, >& "dev$null" )

          print (xphot, yphot, > "cord.out")
          radd = min( (ix2-ix1)/2, (iy2-iy1)/2 )

          phot (cutout,"", , interactive=no, update=yes, coords="cord.out", 
          output="f", plotfile="", 
          scale=1., fwhmpsf=2.5, radplots=no, verify=no, 
          emission=yes, sigma=INDEF, datamin=INDEF, 
          datamax=INDEF, noise="poisson", ccdread="", gain="", readnoise=0., 
          epadu=1., exposure="", airmass="", filter="", obstime="", itime=1., 
          xairmass=INDEF, ifilter="INDEF", otime="INDEF", calgorithm="none",
          cbox=5., cthreshold=0., minsnratio=1., cmaxiter=10, maxshift=1., 
          clean=no, rclean=1., rclip=2., kclean=3., mkcenter=no, 
          salgorithm="median", annulus=10., dannulus=10., skyvalue=0., 
          smaxiter=10, sloclip=0., shiclip=0., snreject=50, sloreject=3., 
          shireject=3., khist=3., binsize=0.1, smooth=no, rgrow=0., 
          mksky=no, weighting="constant", apertures=radd, zmag=25., 
          mkapert=no,  
          verbose=no, graphics="stdgraph", display="stdimage", icommands="",
          gcommands="")
       
          phot.interactive = no
          phot.verify = no

          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
        mylist2 = "fff"

	while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, btakeout, npix )
}
           btakeout = btakeout / npix

#           print (npix)

# subtract BKGND mean from cutout -> ssi

          imarith ( cutout, "-", btakeout, ssi, title="", divzero=0., hparams="", pixtype="real", calctype="real", verbose-, noact- )

# compute initial rotation center

	  bgx = (x2 - x1) / 2 + dc
	  bgy = (y2 - y1) / 2 + dc

# iterate on minimization until convergence

	  cflag = no
	  bcount = 0

	  while ( !cflag ) {

	    bcount += 1

	    xcen = bgx
	    ycen = bgy

# spiral walk rotation center through 3x3 box

	    for ( kk=1; kk<=9; kk+=1 ) {

	      if ( kk == 2 || kk == 8 || kk == 9 )  xcen = xcen-search
	      if ( kk == 4 || kk == 5 )  xcen = xcen+search
      	      if ( kk == 3 )  ycen = ycen+search
	      if ( kk == 6 || kk == 7 )  ycen = ycen-search

# Rotate bkgnd-subtracted image, subtract original, and take abs value of
# difference: ssi -> rot -> sub -> sqr

	      rotate  ( ssi, rot, 180., xin=xcen, yin=ycen, xout=xcen,
                        yout=ycen, ncols=INDEF, nlines=INDEF, interpolant="linear",
                        boundary="nearest", constant=0., nxblock=256, nyblock=256,
		        verbose- )
              imarith ( rot, "-", ssi, sub, title="", divzero=0.,
		        hparams="", pixtype="real", calctype="real", verbose-,
		        noact- )
	      imfunction ( sub, sqr, "abs", verbose- )

# get mean value of abs. val. of resid. image: takeout

           xphot = xcen
	   yphot = ycen

           delete ( "cntr", verify-, >& "dev$null" )
           delete ( "cord.out", verify-, >& "dev$null" )
           delete ( "f", verify-, >& "dev$null" )
           delete ( "ff.tab", verify-, >& "dev$null" )
           delete ( "fff", verify-, >& "dev$null" )

          print (xphot, yphot, > "cord.out")
          radd = min( (ix2-ix1)/2, (iy2-iy1)/2 )

        if ( radd+xphot >= x2-x1+2*dc  ) {radd = x2-x1+2*dc - xphot}
        if ( radd+yphot >= y2-y1+2*dc  ) {radd = y2-y1+2*dc - yphot}
        if ( xphot - radd <= 1  ) {radd = xphot - dc}
        if ( yphot - radd <= 1  ) {radd = yphot - dc}

          phot (sqr, "", coords="cord.out", output="f", plotfile="", 
          scale=1., apertures=radd, calgorithm="none", zmag=25.)

          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
           mylist2 = "fff"
	   while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, takeout, npix )
}

           takeout = takeout / npix


# Takeout is the flux from the rotated and subtracted image	  

# Do the minimization: update best asymm value on each pass.

  	      if ( kk == 1 ) {
	        symm = takeout
	        centr = symm

                imstat ( sqr//section, fields="stddev", lower=INDEF,
		         upper=INDEF, binwidth=0.1, format-, >> tmpdat )

                mylist2 = tmpdat
                jj = fscan ( mylist2, stdd )
     	        delete ( tmpdat, verify-, >& "dev$null" )
	      }

	      tsymm = takeout

	      if ( tsymm < symm ) {
	        symm = tsymm
	        bgx = xcen
	        bgy = ycen
	      }
# clean up
              imdelete ( rot, verify-, >& "dev$null" )
              imdelete ( sub, verify-, >& "dev$null" )
              imdelete ( sqr, verify-, >& "dev$null" )

	    }					# close for kk=1,9 spiral walk

	    if ( symm == centr ) cflag = yes

	  }					# close for while (!cflag)

	  takeout = symm
          imdelete ( cutout, verify-, >& "dev$null" )
          imdelete ( ssi, verify-, >& "dev$null" )
#
#  END BACKGROUND

	flprc; flprc

	  ix1 = int ( x - 250 )
          ix2 = int ( x + 249 )
	  iy1 = int ( y - 250 )
	  iy2 = int ( y + 249 )

	  x = 251
	  y = 251
	  rad = 251

	  if (xo-250 < 1 || xo+250>xdim-1 || yo-250 < 1 || yo+250 > ydim-1) {
                     nrad = min(xo,yo,xdim-xo,ydim-yo)
                     nrad = nrad - 1
	              ix1 = int ( xo - nrad )
                      ix2 = int ( xo + nrad-1 )
	              iy1 = int ( yo - nrad )
	              iy2 = int ( yo + nrad-1)
		      x = nrad+1; y = nrad+1
                      rad = nrad

		      size = 2*nrad}
#		      print (size)

	  section = "["//ix1//":"//ix2//","//iy1//":"//iy2//"]"
#	print (btakeout)
        imarith (image, "-", btakeout , "test", title="", divzero=0., hparams="", pixtype="real", calctype="real", verbose=no, noact=no) 

             	imcopy ("test"//section, "im",  verbose- )

#   Make a (small) copy of the image to work with, called im

	kk = 1

	cflag = no
	iflag = no
	oflag = no
        halfcon = no

	tsb = 0

        imcntr ("im", x, y, cboxsize=15, > "cntr")

        mylist2 = "cntr"

	while ( fscan (mylist2, junk1, junk2, x, junk3, y) !=EOF )  
           delete ( "cntr", verify-, >& "dev$null" )


#   FIND THE CENTRAL SURFACE BRIGHTNESS

# 	BEGIN THE NEW SEARCH FOR THE RADII WHERE THE SURFACE BRIGHTENSS IS
# 	ONE HALF OF THE MAX.

	while ( !cflag ) {

#	print (kk)

	ix1 = int ( x - kk )
	ix2 = int ( x + kk )
	iy1 = int ( y - kk )
	iy2 = int ( y + kk )


        if ( x+kk+1 >= size ) {cflag =yes
			      kk = size - x - 1
			      flag = 1}

        if ( y+kk+1 >= size ) {cflag =yes
			      kk = size - y -1
			      flag = 1}

        if ( x-kk-1 < 0 ) {cflag =yes
			   kk = x - 1
			   flag = 1}

        if ( y-kk-1 < 0 ) {cflag =yes
	  		   kk = y - 1
			   flag = 1}
       
        if (kk > 249) {cflag = yes}

        section = "["//ix1//":"//ix2//","//iy1//":"//iy2//"]"

#	print (kk, flag, size, x, y)

#PHOT START

	    xphot = x
	    yphot = y

#	GET THE OUTER SB

           delete ( "cord.out", verify-, >& "dev$null" )
           delete ( "f", verify-, >& "dev$null" )
           delete ( "ff.tab", verify-, >& "dev$null" )
           delete ( "fff", verify-, >& "dev$null" )

          print (xphot, yphot, > "cord.out")

          radd = kk

          phot ("im", "", coords="cord.out", output="f", plotfile="", 
          scale=1., apertures=radd, calgorithm="none", zmag=25.)

          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
         mylist2 = "fff"
	 while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, osb, npix )
}
#           print (osb, npix)
#           b = 1/npix
#OSB - the outer part of the sb profile starting at pixel 2
#           print (osb)


# GET THE FLUX FROM THE INNER REGION.

	section = "["//ix1+1//":"//ix2-1//","//iy1+1//":"//iy2-1//"]"

           delete ( "cord.out", verify-, >& "dev$null" )
           delete ( "f", verify-, >& "dev$null" )
           delete ( "ff.tab", verify-, >& "dev$null" )
           delete ( "fff", verify-, >& "dev$null" )

          print (xphot, yphot, > "cord.out")

          radd = kk-1
	  
          phot ("im", "", coords="cord.out", output="f", plotfile="", 
          scale=1., apertures=radd, calgorithm="none", zmag=25.)

          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
         mylist2 = "fff"
	 while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, isb, nipix )
}
#         print (isb, npix)

# isb - calculcate the inner part (kk-1) of the profile, starting at 1

           otsb = tsb

#	   print (npix, nipix, osb, npix)

         if (npix < 1) {npix = 1; osb = 1; flag = 1; cflag = yes}

	tsb = ((osb - isb) / (npix - nipix)) / (osb/npix)

#       Remove lower # to print the Petrosian radii (tsb)

#	print ( kk, tsb, osb, tsb * npix / osb  )

#       Find the total radius, by the Petrosian radii at eta = 0.2

        t = tsb / (osb/npix)

	if ( tsb < 0.2)  {
       
        cflag = yes

#INTERPOLATE THE 0.2 PETROSIAN RADIUS

        kk = kk - (1/(otsb - tsb))*(0.2 - tsb) 
}

#        print (tsb, kk, cflag)
   
	kk = kk + 1

	flprc
	flprc
}
 
#End the CFLAG for computing the eta = 0.2 petrosian radius

	cflag = no	

	kk = kk - 1

	if (kk < 2) {kk = 1; flag = 1} 

#       Bigrad = Total radii, dafault = 1.5 * eta = 0.2 Petrosian.
	bigrad =  1.5 * kk

#	FINISHED WITH FINDING THE RADIUS, NOW FIND THE CONCENTRATIONS.
#	FIND TOTAL FLUX WITHIN 1.5 * eta radii

	kk = bigrad

	ix1 = int ( x - kk )
	ix2 = int ( x + kk )
	iy1 = int ( y - kk )
	iy2 = int ( y + kk )

        if ( x+kk+1 > size ) {cflag =yes
			  kk = size - x - 2 - dc
                          bigrad = size - x - 2 - dc
			  flag = 1}

        if ( y+kk+1 > size ) {cflag =yes
			   kk = min(kk, size - y -2 - dc)
                           bigrad = min(bigrad, size - y - 2 - dc)
			   flag = 1}

        if ( x-kk-1 < 0 ) {cflag =yes
			  kk = min(kk, x - 2 - dc)
                          bigrad = min(bigrad, x - 2 - dc)
			  flag = 1}

        if ( y-kk-1 < 0 ) {cflag =yes
			   kk = min(kk, y - 2 - dc)
                           bigrad = min(bigrad, y - 2 - dc)
			   flag = 1}


	section = "["//ix1//":"//ix2//","//iy1//":"//iy2//"]"

           delete ( "cord.out", verify-, >& "dev$null" )
           delete ( "f", verify-, >& "dev$null" )
           delete ( "ff.tab", verify-, >& "dev$null" )
           delete ( "fff", verify-, >& "dev$null" )

          print (xphot, yphot, > "cord.out")
          radd = kk

          phot ("im", "", coords="cord.out", output="f", plotfile="", 
          scale=1., apertures=radd, calgorithm="none", zmag=25.)

          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
         mylist2 = "fff"
	 while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, flux, npix )
}

#	print (flux, radd)

	if (flux < -1) {print (flux, npix, xphot, yphot, radd)
       }
	

#flux is the total flux of the galaxy within the eta = 0.2 radius

#	Begin computing concentrations

	ll = 2

	tflux = 0

	while ( !iflag ) {

	ix1 = int ( x - ll )
	ix2 = int ( x + ll )
	iy1 = int ( y - ll )
	iy2 = int ( y + ll )

        if ( x+ll+1 > size ) {cflag =yes
			  ll = size - x - 1
			  flag = 1}

        if ( y+ll+1 > size ) {cflag =yes
			   ll = size - y -1
			   flag = 1}

        if ( x-ll-1 < 0 ) {cflag =yes
			   ll = x - 1
			   flag = 1}

        if ( y-ll-1 < 0 ) {cflag =yes
			   ll = y - 1
			   flag = 1}
 
        if (ll > 249) {cflag = yes
                       irad = 249
                       iflag = yes
                       flag = 1}

	section = "["//ix1//":"//ix2//","//iy1//":"//iy2//"]"

	otflux = tflux

	xphot = x
	yphot = y

           delete ( "cord.out", verify-, >& "dev$null" )
           delete ( "f", verify-, >& "dev$null" )
           delete ( "ff.tab", verify-, >& "dev$null" )
           delete ( "fff", verify-, >& "dev$null" )

          print (xphot, yphot, > "cord.out")
          radd = ll

          phot ("im", "", coords="cord.out", output="f", plotfile="", 
          scale=1., apertures=radd, calgorithm="none", zmag=25.)

          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
         mylist2 = "fff"
	 while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, tflux, npix )
}

     	if (tflux < -1) {iflag = yes
                         irad = 10}

#       Determines the inner radii for concentration, default = 0.2

	if ( tflux > 0.2 * flux) {
	iflag = yes

#	print (otflux, tflux, flux)

        irad = ll - ((0.2*flux - tflux)/(otflux - tflux)) 
	otflux = 0}

	if (ll == 0) {
	iflag = yes}

	ll = ll + 1
}
#End finding the inner 0.2 radius

        tflux = 0
	iflag = no

#	FIND THE OUTER RADIUS

	ll = irad - 1

	while ( !oflag ) {

	ix1 = int ( x - ll )
	ix2 = int ( x + ll )
	iy1 = int ( y - ll )
	iy2 = int ( y + ll )

        if ( x+ll+2 > size ) {cflag =yes
			    ll = size - x - 1
			    flag = 1}

        if ( y+ll+2 > size ) {cflag =yes
			    ll = size - y -1
			    flag = 1}

        if ( x-ll-1 < 0 ) {cflag =yes
			   ll = x - 1
			   flag = 1}

        if ( y-ll-1 < 0 ) {cflag =yes
			   ll = y - 1
			   flag = 1}

        if (ll > 249) {cflag = yes
                       oflag = yes
                       orad = 248}

	section = "["//ix1//":"//ix2//","//iy1//":"//iy2//"]"

	xphot = x
	yphot = y

           delete ( "cord.out", verify-, >& "dev$null" )
           delete ( "f", verify-, >& "dev$null" )
           delete ( "ff.tab", verify-, >& "dev$null" )
           delete ( "fff", verify-, >& "dev$null" )

          print (xphot, yphot, > "cord.out")
          radd = ll

          phot ("im", "", coords="cord.out", output="f", plotfile="", 
          scale=1., apertures=radd, calgorithm="none", zmag=25.)

          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
         mylist2 = "fff"
	 while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, tflux, npix )
}
	if (tflux < -1) {oflag = yes
	                 orad = 10}
	

#       Determine half-light radius

        if (!halfcon) {

        if ( tflux > 0.5 * flux) {
          halfcon = yes
          halfr = ll - ((0.5*flux - tflux)/(otflux - tflux))}}

# Deteremines the outer radii for concentration, default = 0.8
   
	if ( tflux > 0.8 * flux) {
	oflag = yes
	orad = ll
        orad = ll - ((0.8*flux - tflux)/(otflux - tflux))}

	if ( ll == 0 ) {
	oflag = yes}

	ll = ll + 1 } # end oflag

        flprc
        flprc

	oflag = no 

	if (orad < 2 || irad < 1) {orad = 2
                                  irad = 2}      

	econc = 5*log10 ( (orad) / (irad) ) - 5*log10 ( (orad-0.5) / (irad+0.5) )
	conc = 5 * log10 (orad / irad)
        conc = real(int(1000. * conc))/1000.
        econc = real(int(1000. * econc))/1000.
        orad = real(int(1000. * orad))/1000.
        irad = real(int(1000. * irad))/1000.
        bigrad = real(int(1000. * bigrad))/1000.
        halfr = real(int(1000. * halfr))/1000.


	flprc
	flprc
	flprc

#	print ('begin asymmetry')

#############################################################################
#                                                                           #
#
#                      Begin Asymmetry Source computation
#                        Note 'takeout' computed earlier                    #
#                                                                           #
#############################################################################
          
          rad = bigrad
	  rad0 = rad

	  if ( x+rad+dc+1>size || x-rad-dc-1<0 || y+rad+dc+1>size || y-rad-dc-1<0 ) {
#	     if ( x+rad0+dc+1>size ) newrad = min( rad-dc-1, x-2)
	newrad = size

             if ( x+rad0+dc+1>size ) newrad = size - x - dc - 2
	     if ( y+rad0+dc+1>size ) newrad = min(newrad, size - y - dc - 2)
	     if ( x-rad0-dc-1<0 ) newrad = min(newrad, x - dc - 2)
	     if ( y-rad0-dc-1<0 ) newrad = min(newrad, y - dc - 2)

#	     if ( y+rad0+dc+1>size ) newrad = min( rad-dc-1, y-2)
# 	     if ( x-rad0-dc-1<0 )    newrad = min( rad-dc-1, x-2)
#	     if ( y-rad0-dc-1<0 )    newrad = min( rad-dc-1, y-2)
        	rad = newrad
	     print ( 'WARNING: user radius is off image; ' )
	     print ( 'radius reset from ',rad0,' to ',rad )

	  }

# define extraction box and extract

	  xorg = x
	  yorg = y

	  ixx1 = int ( xorg-rad)
          ixx2 = int ( xorg+rad)

          iyy1 = int ( yorg-rad)
          iyy2 = int ( yorg+rad)

	  psect = "["//ixx1-dc+1//":"//ixx2+dc-1//","//iyy1-dc+1//":"//iyy2+dc-1//"]"
	  section = "["//dc//":"//ixx2-ixx1+dc//","//dc//":"//iyy2-iyy1+dc//"]"
             imcopy ("im"//psect, ssi,  verbose- )
# subtract BKGND

# create abs value image and get statistics for denominator value

	  imfunction ( ssi, imsqr, "abs", verbose- )

	  xphot = dc + rad 
	  yphot = dc + rad

           delete ( "cord.out", verify-, >& "dev$null" )
           delete ( "f", verify-, >& "dev$null" )
           delete ( "ff.tab", verify-, >& "dev$null" )
           delete ( "fff", verify-, >& "dev$null" )

          print (xphot, yphot, > "cord.out")

          radd = rad

        if ( xphot+radd+1 > 2*rad + 2*dc ) {cflag =yes
			    radd = 2*rad + 2*dc - xphot
			    flag = flag + 1}

        if ( yphot+radd+1 > 2*rad + 2*dc ) {cflag =yes
			    radd = 2*rad + 2*dc - yphot
			    flag = flag + 1}
			
        if ( xphot-radd-1 < 0 ) {cflag =yes
			   radd = rad + dc - 1
			   flag = flag + 1}

        if ( yphot-radd-1 < 0 ) {cflag =yes
			   radd = rad + dc - 1
			   flag = flag + 1}

          phot (imsqr, "", coords="cord.out", output="f", plotfile="", 
          scale=1., apertures=radd, calgorithm="none", zmag=25.)

          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
         mylist2 = "fff"
	 while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, bottom, npix )
}
           bottom = bottom / npix

#print (npix)
	  delete ( tmpdat, verify-, >& "dev$null" )

          imstat ( imsqr//section, fields="stddev", lower=INDEF, upper=INDEF,
        	   binwidth=0.1, format-, >> tmpdat )

	  mylist2 = tmpdat

	  jj = fscan ( mylist2, istd )

	  delete ( tmpdat, verify-, >& "dev$null" )
	  imdelete ( imsqr, verify-, >& "dev$null" )

# compute initial rotation center

	  x = rad + dc+1
	  y = rad + dc+1

# iterate on minimization until convergence

	  cflag = no
	  ccount = 0

	  while ( !cflag ) {

	    ccount += 1

	    xcen = x
	    ycen = y

# spiral walk rotation center through 3x3 box

	    for ( kk=1; kk<=9; kk+=1 ) {

	      if ( kk == 2 ) xcen = xcen-search
	      if ( kk == 8 ) xcen = xcen-search
	      if ( kk == 9 ) xcen = xcen-search
	      if ( kk == 4 || kk == 5 )  xcen = xcen+search
      	      if ( kk == 3 )  ycen = ycen+search
	      if ( kk == 6 || kk == 7 )  ycen = ycen-search

# Rotate bkgnd-subtracted image, subtract original, and take abs value of
# difference: ssi -> rot -> sub -> sqr

              rotate ( ssi, rot, 180., xin=xcen, yin=ycen, xout=xcen,
		       yout=ycen, ncols=INDEF, nlines=INDEF, interpolant="linear", boundary="nearest", constant=0., nxblock=256, nyblock=256, verbose- )

              imarith ( rot, "-", ssi, sub, title="", divzero=0.,
		        hparams="", pixtype="real", calctype="real", verbose-, noact- )
 	      imfunction ( sub, sqr, "abs", verbose- )


# get mean value of abs. val. of resid. image: top

           delete ( "cord.out", verify-, >& "dev$null" )
           delete ( "f", verify-, >& "dev$null" )
           delete ( "ff.tab", verify-, >& "dev$null" )
           delete ( "fff", verify-, >& "dev$null" )
	   
#	   print (xcen, ycen, radd, rad, dc)

        if ( xcen+radd+1 > 2*rad + 2*dc ) {cflag =yes
			    radd = 2*rad + 2*dc - xcen -1
			    flag = flag + 1}

        if ( ycen+radd+1 > 2*rad + 2*dc ) {cflag =yes
			    radd = 2*rad + 2*dc - ycen -1
			    flag = flag + 1}
			
        if ( xcen-radd-1 < 0 ) {cflag =yes
			   radd = xcen - 1
			   flag = flag + 1}

        if ( ycen-radd-1 < 0 ) {cflag =yes
			   radd = ycen - 1
			   flag = flag + 1}
          xphot = xcen
 	  yphot = ycen

          print (xphot, yphot, > "cord.out")

          phot (sqr, "", coords="cord.out", output="f", plotfile="", 
          scale=1., apertures=radd, calgorithm="none", zmag=25.)

          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
         mylist2 = "fff"
	 while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, top, npix )
}
#	print (npix)
           top = top / npix

# Do the minimization: update best asymm value on each pass.

	      if ( kk == 1 ) {
	        symm = (top - takeout) / (bottom)
	        centr = symm
	      }

	      tsymm = (top - takeout) / (bottom)

	      if ( tsymm < symm ) {
	        symm = tsymm
	        x = xcen
	        y = ycen
	      }

# Error computation
	      error = (stdd / (2*(bottom+istd)))

# clean up
              imdelete ( rot, verify-, >& "dev$null" )
              imdelete ( sub, verify-, >& "dev$null" )
              imdelete ( sqr, verify-, >& "dev$null" )

	    }				# close for kk=1,9 spiral walk

	    if ( symm == centr ) cflag = yes

	  }    # closing for while (!cflag)

          imdelete ( cutout, verify-, >& "dev$null" )
          imdelete ( ssi, verify-, >& "dev$null" )

# update center

	  xp = xcen + ixx1 - dc
	  yp = ycen + iyy1 - dc

#print (xcen, ixx1, dc)
#print (xp, yp, xorg, yorg)

# Reformat symm and error to take 3 sig. dig., and output

	  symm = real(int(1000. * symm))/1000.
	  error = real(int(1000. * error))/1000
	  xorg = real(int(1000. * xorg))/1000
	  yorg = real(int(1000. * yorg))/1000
	  xp = real(int(1000. * xp))/1000
          yp = real(int(1000. * yp))/1000
	  rad = real(int(1000. * rad))/1000

	  x = xp
	  y = yp

 #        print ('Begin Computation of S')

	flprc
	flprc
	flprc

#########################################################################
#
#
#      Start the high spatial freq. filtering
#
#        
#########################################################################

	imdelete ( cutout, verify-, >& "dev$null" )
	imdelete ( ssi, verify-, >& "dev$null" )
	imdelete ( rot, verify-, >& "dev$null" )
	imdelete ( sub, verify-, >& "dev$null" )
	imdelete ( sqr, verify-, >& "dev$null" )
	imdelete ( imsqr, verify-, >& "dev$null" )
	delete ( tmpdat, verify-, >& "dev$null" )

	if (radd < 1) {radd = 3}

#CHANGE THIS  RAD/ 6 TO RAD / ?
	        sm = radd / 6

# following 3 lines added by jeff b to fix smoothing problem 20020829

		sm = int (sm)
		jgv = sm / 2
		if ( jgv == int(jgv)) { sm = sm + 1}	#changes even numbers to odds


############################################################################
# Source calculation
############################################################################

          imcopy ("im", ssi,  verbose- )

#          imarith (cutout, "-", btakeout, ssi, title="",  divzero=0., 
#        	    hparams="", pixtype="", calctype="",verbose-, noact- )
# subtract BKGND
 
	  boxcar (ssi, rot, sm, sm, boundary="nearest", constant=0.)

	  imarith (ssi,"-", rot, "c1", title="", divzero=0., hparams="",
          pixtype="", calctype="", verbose=no, noact=no)

	  imreplace ("c1", 0., imaginary=0., lower=-50000., upper=0., radius=0.)
	  ix1 = int ( x - radd/20)
          ix2 = int ( x + radd/20)
          iy1 = int ( y - radd/20)
          iy2 = int ( y + radd/20)

	  psect = "["//ix1//":"//ix2//","//iy1//":"//iy2//"]"

	  imreplace ("c1"//psect, 0., imaginary=0., lower=0., upper=99999., radius=0.)

#Flux from the subtracted image

           delete ( "cord.out", verify-, >& "dev$null" )
           delete ( "f", verify-, >& "dev$null" )
           delete ( "ff.tab", verify-, >& "dev$null" )
           delete ( "fff", verify-, >& "dev$null" )
	   
	   xphot = x
	   yphot = y

          print (xphot, yphot, > "cord.out")

        if ( x+radd+1 > size ) {cflag =yes
			      radd = radd - x - 1
			      flag = flag + 1}

        if ( y+radd+1 > size ) {cflag =yes
			      radd = radd - y -1
			      flag = flag + 1}

        if ( x-radd-1 < 0 ) {cflag =yes
			      radd = x - 1
			      flag = flag + 1}

        if ( y-radd-1 < 0 ) {cflag =yes
			      radd = y -1
			      flag = flag + 1}

          phot ("c1", "", coords="cord.out", output="f", plotfile="", 
          scale=1., apertures=radd, calgorithm="none", zmag=25.)

          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
         mylist2 = "fff"
	 while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, top, npix )
}

           top = top / npix

#Flux from the total image

           delete ( "cord.out", verify-, >& "dev$null" )
           delete ( "f", verify-, >& "dev$null" )
           delete ( "ff.tab", verify-, >& "dev$null" )
           delete ( "fff", verify-, >& "dev$null" )

          print (xphot, yphot, > "cord.out")

          phot (ssi, "", coords="cord.out", output="f", plotfile="", 
          scale=1., apertures=radd, calgorithm="none", zmag=25.)

          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
         mylist2 = "fff"
	 while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, bottom, npix )
}
           bottom = bottom / npix

          	imdelete ( cutout, verify-, >& "dev$null" )
          	imdelete ( ssi, verify-, >& "dev$null" )
          	imdelete ( rot, verify-, >& "dev$null" )
          	imdelete ( "c1", verify-, >& "dev$null" )
          	imdelete ( "t1", verify-, >& "dev$null" )
          	imdelete ( "c2", verify-, >& "dev$null" )

# get BKGND mean (btakeout) from the medianed image.
#Take out the background contribution to the S parameter

	bsize = min ( (x2-x1), (y2-y1) )

#new stuff below for the faster method.

		ix1 = int ( x1 )
		ix2 = int ( x2 )
		iy1 = int ( y1 )
		iy2 = int ( y2 )

		bsect = "["//ix1-1//":"//ix2+1//","//iy1-1//":"//iy2+1//"]"

          	imcopy (image//bsec, "cutout",  verbose- )
          	imarith ("cutout", "-", btakeout, ssi, title="",  divzero=0., hparams="", pixtype="", calctype="",verbose-, noact- )

	  	flprc
	  	flprc

	  	boxcar (ssi, rot, sm, sm, boundary="nearest", constant=0.)

	  	imarith (ssi,"-", rot, "c1", title="", divzero=0., hparams="", pixtype="", calctype="", verbose=no, noact=no)

	  	imreplace ("c1", 0., imaginary=0., lower=-50000., upper=0., radius=0.)

           delete ( "cord.out", verify-, >& "dev$null" )
           delete ( "f", verify-, >& "dev$null" )
           delete ( "ff.tab", verify-, >& "dev$null" )
           delete ( "fff", verify-, >& "dev$null" )

	   xphot =  (x2 - x1)/2
	   yphot =  (y2 - y1)/2 

          print (xphot, yphot, > "cord.out")

          radd = min( (x2-x1)/2, (y2-y1)/2 )
	  radd = radd - 1

        if ( xphot+radd+1 > bsize ) {cflag =yes
			      radd = radd - xphot - 1
			      }

        if ( yphot+radd+1 > bsize ) {cflag =yes
			      radd = radd - yphot -1
			      }

        if ( xphot-radd-1 < 0 ) {cflag =yes
			      radd = xphot - 1
			      }

        if ( yphot-radd-1 < 0 ) {cflag =yes
			      radd = yphot - 1
			      }

          phot ("c1", "", coords="cord.out", output="f", plotfile="", 
          scale=1., apertures=radd, calgorithm="none", zmag=25.)

          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
         mylist2 = "fff"
	 while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, ttakeout, npix )
}
           ttakeout = ttakeout / npix

#Delete files

          imdelete ( "cutout", verify-, >& "dev$null" )
          imdelete ( cutout, verify-, >& "dev$null" )
          imdelete ( ssi, verify-, >& "dev$null" )
          imdelete ( rot, verify-, >& "dev$null" )
          imdelete ( "c1", verify-, >& "dev$null" )
          imdelete ( "t1", verify-, >& "dev$null" )
          imdelete ( "c2", verify-, >& "dev$null" )

# Reformat s and error to take 3 sig. dig., and output

          s = (top - ttakeout) / bottom 
	  se = 10*(sqrt(abs(npix*(top-ttakeout))))/(npix*(bottom+istd))
	  s = 10 * (real(int(1000. * s))/1000.)	  
	  se = real(int(1000. * se))/1000.

#END the computation
 
	delete ( tmpdat, verify-, >& "dev$null" )

# Reformat symm and error to take 3 sig. dig., and output

	  symm = real(int(1000. * symm))/1000.
	  error = real(int(1000. * error))/1000
	  xorg = real(int(1000. * xorg))/1000
	  yorg = real(int(1000. * yorg))/1000
	  xp = real(int(1000. * (xo+xp - xorg - 1)))/1000
          yp = real(int(1000. * (yo+yp - yorg - 1)))/1000 
	  rad = real(int(1000. * rad))/1000
	  econc = real(int(1000. * econc))/1000

#renormalize the flux 
	  flux = flux / 12000
	  be = bigrad
	  be = real(int(1000. * be))/1000.
	  flux = real(int(1000. * flux))/1000.

   	print ( image, " ", xp, yp, conc, econc, symm, error, s, se, bigrad, halfr, flux, flag, >> ofile)

	print ( image, " ", conc, econc, symm, error, s, se, bigrad, halfr, flux, flag )

	flag = 0

        flprc
        flprc

}
 
end
 

