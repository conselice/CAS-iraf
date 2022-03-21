## Version of CAS running on a list of objects located on different tiles
#  Programs (mtcaspkg.cl,  mtcaswrap.cl, mtcas.cl)
#
# 1) The program mtcaspkg.cl is  a wrapper package where task
#   mtcaswrap.cl runs the  CAS code mtcas.cl (by C  Conselice) on an
#   input list of objects found on multiple tiles. Format of input file is 
#           (object id, image,segmap, x,y, x1,x2,y1,y2)
#           where x y =center,  and x1,x2,y1,y2 defines sky box
# 2) To run  
# a) Edit  login.cl.
#   set imtype =fits 
#   set     scrpath = "/data/azur8/irafsc/"
#   task  $mtcaspkg = scrpath$mtcaspkg.cl 
#
# b) Start iraf. 
#    cl> mtcaspkg    ; load package
#    cl> cd /data/gems   (= where batch-file1 and tiles are)
#    cl> mtcaswrap.cl  inputfile outputfile
#   where inputfile =(object id, image,segmap, x,y, x1,x2,y1,y2)
#     
# 
######################################################################
# CAS Pg version -1 by C.Conselice to compute radii, asymmetry and 
# concentration indexes, fluxes.  #  
#  Contact conselice@gmail.com for questions
#
# Change x1,x2,y1,y2 to be the coordinates of the background light
#
# x-center,y-center are the centers of the galaxies under study, background
# is the background level of the image.
#
# Note the following conditions on input values:
#       x1 - del > 0
#       y1 - del > 0
#       x2 + del < xdim
#       y2 + del < ydim
# If the input radius is plus coordinates exceed image size then 
# is rest to maximum size that can fit on image.
#  
#
###SJ/JDbegin
#multi-object cas, oct 16, 2004
#procedure cas ( imgfile, outfile )

procedure mtcas (infile, outfile)

string infile = ""{prompt="Multiple Object CAS -- Input file: id, image,  wmap_image, segmap, x,y" } 

string outfile = ""{prompt="Outfile File Name" }


###
#string imgfile = ""	{ prompt = "Input file: image x y background size" }
#string outfile = ""	{ prompt = "Output file: image, x, y, conc, orad, irad, petr, halfr, flag" }
###SJ/JDend

real search 	{ prompt = "Pixel step-size for minimization search" }
int dc 		{ prompt = "Maximum expected center wander for minimization" }

bool vbo1 = no	{ prompt = "Level 1 diag: obj. and calc. step" }
bool vbo2 = no	{ prompt = "Level 2 diag: minimization iterations" }

struct *mylist1		{ prompt = "wiz" }
struct *mylist2		{ prompt = "bang" }
struct *mylist3         { prompt = "zoom" }

begin



	int ix1, ix2, iy1, iy2
	int ii, jj, flag, red, xd, yd, xseg, yseg
	string image, timage, wimage, imaget, sub, st, junk1, junk2, junk3
	string imsqr, segsec
        string segmap
	string sqr, ifile, ofile
	string rot, section, oname
	string z, simage#gini stuff below
#gini stuff
        real numc,aveflux,gini1,gini2,gini

	real xt, yt, inten, dis, fdis, tdis, tabflux, c20, xm,  ym
        real bigrad, radius, ll, se, xo, xc, yc, yo, bsize
        real sflag, segval
	real  t, bg, csb, isb, osb, size, otflux
	real flux, kk, cc, id, otsb, econc, wobj
	real npix, tsb, nipix, tflux, halfr, si
	real irad, orad, conc, name
	string  cutout, ssi, tmpdat
	real symm, top, bottom, tsymm,gsymm, centr, xphot, yphot, radd
	real bgx, bgy, xp, yp, rad0, e, be, jgv
	real xcen, ycen, x1, x2, y1, y2, rad, xorg, yorg, stddev
	real wtakeout, otakeout, btakeout, takeout, stdd, istd, error
	bool cflag, iflag, oflag, halfcon, dflag
	int ixx1,ixx2,iyy1,iyy2
	int xdim,ydim,ccount,bcount
	string psect

	int dd
	string bsect
	real sm, s
        real eta8, eta5, eta2, grw2, grw5, grw8, eta2t2, grw8t2, ttakeout
noao
digiphot
apphot
tables
ttools
artdata

###SJ/JDbegin
# Check if necessary packages are loaded
if (!deftask("apphot")) {
        print ("ERROR: Load noao.digiphot.apphot packages.")
        bye
}
if (!deftask("artdata")) {
        print ("ERROR: Load artdata package.")
        bye
}
if (!deftask("ttools")) {
        print ("ERROR: Load tables.ttools packages.")
        bye
}
# Load up packages
#        noao
#        digiphot
#        apphot
#        artdata 
#        tables
#        ttools
###SJ/JDend

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

### BEGIN  EDITS
#initialize inputs
	ifile = infile
	list = ifile
	ofile = outfile
#	wimage = wile
#	ifile = imgfile
#	mylist1 = ifile

        print ( '#  Name  x        y      C    E(C)    A    E(A)   S   E(S)  M_20 Gini Petr   HalfR  Flux   inflag casflag', >> ofile )

 	while (fscan(list, timage, id, x, y, wimage, simage, sflag) !=EOF) {

                flprc
                flprc    

	  imget ( simage, "i_naxis1" )
	  xdim = int(imget.value)
	  imget ( simage, "i_naxis2" )
	  ydim = int(imget.value)

       cflag = no

       segmap = simage

#calculate the new background based on the x,y guess


        x1 = x 
        x2 = x 
        y1 = y 
        y2 = y

      cc = 0

	  while ( !cflag ) {

        cc = cc + 1

	    for ( kk=1; kk<=9; kk+=1 ) {

	      if ( kk == 2 || kk == 8 || kk == 9 )  x1 = x1-cc*1
	      if ( kk == 4 || kk == 5 )  x1 = x1+cc*1
      	      if ( kk == 3 )  y1 = y1+1
	      if ( kk == 6 || kk == 7 )  y1 = y1-1

             x2 = x1 + 10
             y2 = y1 + 10

          if (x1 <= 1) {x1 = x + 75; x2 = x1 - 10}
          if (y1 <= 1) {y1 = y + 75; y2 = y1 - 10}

          if (x2 >= xdim) {x1 = xdim - 10; x2 = xdim - 1}
          if (y2 >= ydim) {y1 = ydim - 10; y2 = ydim - 1}

 
        ix1 = int (x1)
        ix2 = int (x2)
        iy1 = int (y1)
        iy2 = int (y2)

# print (kk, ix1, ix2, iy1, iy2)


	section = "["//ix1//":"//ix2//","//iy1//":"//iy2//"]"
#print (xdim, ydim, section)
                imstat ( segmap//section, fields="mean", lower=INDEF,
		         upper=INDEF, binwidth=0.1, format-, >> 'tmpdat' )

                mylist2 = 'tmpdat'
                jj = fscan ( mylist2, segval )
     	        delete ('tmpdat', verify-, >& "dev$null" )
          
              if (cc > 500) {segval = 0}	      

	      if ( segval == 0 ) {
	        cflag = yes
		      }}}

###end the new background area computation.

### END EDITS

	cutout = mktemp ( "CuT" ) // ".fits"
	ssi = mktemp ( "SsI" ) // ".fits"
	rot = mktemp ( "RoT" ) // ".fits"
	sub = mktemp ( "SuB" ) // ".fits"
	sqr = mktemp ( "SqR" ) // ".fits"
	imsqr = mktemp ( "ISqR" ) // ".fits"
	tmpdat = mktemp ( "TdA" )

	imdelete ( "temp1", verify-, >& "dev$null" )
	imdelete ( "temp2", verify-, >& "dev$null" )
	imdelete ( "cutout.fits", verify-, >& "dev$null" )
	imdelete ( cutout, verify-, >& "dev$null" )
	imdelete ( ssi, verify-, >& "dev$null" )
	imdelete ( rot, verify-, >& "dev$null" )
	imdelete ( sub, verify-, >& "dev$null" )
	imdelete ( sqr, verify-, >& "dev$null" )
	imdelete ( imsqr, verify-, >& "dev$null" )
	imdelete ( "c1.fits", verify-, >& "dev$null" )
	imdelete ( "wcutout.fits", verify-, >& "dev$null" )
	imdelete ( "wseg.fits", verify-, >& "dev$null" )

	delete ( tmpdat, verify-, >& "dev$null" )

	  imdelete ( "CuT*fits*", verify-, >& "dev$null" )
	  imdelete ( "RoT*fits*", verify-, >& "dev$null" )
	  imdelete ( "SsI*fits*", verify-, >& "dev$null" )
	  imdelete ( "SuB*fits*", verify-, >& "dev$null" )
	  imdelete ( "SqR*fits*", verify-, >& "dev$null" )
        image = "test"

	flag = 0
        search = 0.5
        dc = 5
        halfr = 0
        tsb = 0
        
        imdelete (image, verify-, >& "dev$null" )
#        imdel ("test.fits", yes, verify=no, default_acti=yes)
#
#        print ( '#Name           x      y       C     E(C)      A      
#E(A)     S     E(S)    Petr   HalfR   Flux   f', >> ofile )

#
#
#	BEGIN LOOKING AT THE OBJECTS

#        mylist1 = ifile
 
#  Shardha note 
#  To put a new tile=tile.fits ,  do a new symbolic link outside of code
#    ln -s /data/oslo3/good02/z/s1z02A_sci.fits tile.fits 
#
###SJ/JDbegin
#   tima    cache tprint

# timage is set to the input tile image name above in the variable
# initialization section

 

    imarith (timage, "*", 12000 , image, title="", divzero=0., hparams="", pixtype="real", calctype="real", verbose=no, noact=no) 

	  imget ( image, "i_naxis1" )
	  xdim = int(imget.value)
	  imget ( image, "i_naxis2" )
	  ydim = int(imget.value)

	  flprc; flprc

# get image dims

          size = xdim

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

###SJ/JDbegin
#          ix1 = int ( x1)
#          ix2 = int ( x2)
#          iy1 = int ( y1)
#          iy2 = int ( y2)
# these values are initialized above
###SJ/JDend

          if (ix1 - dc <= 0) {ix1 = 1 + dc}
          if (ix2 - dc <= 0) {ix2 = 3 + dc}
          if (iy1 - dc <= 0) {iy1 = 1 + dc}
          if (iy2 - dc <= 0) {iy2 = 3 + dc}

##8/4/03 - here is where the program checks to see if things are too
## high -
          if (ix1 + dc >= xdim) {ix1 = xdim - dc - 2}
          if (ix2 + dc >= xdim) {ix2 = xdim - dc - 1}
          if (iy1 + dc >= ydim) {iy1 = ydim - dc - 2}
          if (iy2 + dc >= ydim) {iy2 = ydim - dc - 1}

	  psect = "["//ix1-dc//":"//ix2+dc//","//iy1-dc//":"//iy2+dc//"]"
	  section = "["//dc//":"//ix2-ix1+dc//","//dc//":"//iy2-iy1+dc//"]"
	  bsect = "["//ix1//":"//ix2//","//iy1//":"//iy2//"]"

#print (psect)

          imcopy ( image//psect, cutout, verbose- )
          imcopy ( wimage//psect, "wcutout", verbose- )

# get BKGND mean (btakeout) from "cutout" image section

          xphot = dc + (ix2 - ix1)/2 
          yphot = dc + (iy2 - iy1)/2 

           delete ( "cord.out", verify-, >& "dev$null" )
           delete ( "f", verify-, >& "dev$null" )
           delete ( "ff.tab", verify-, >& "dev$null" )
           delete ( "fff", verify-, >& "dev$null" )

          print (xphot, yphot, > "cord.out")
          radd = min( (ix2-ix1)/2, (iy2-iy1)/2 ) + 4

          phot (cutout,"", coords="cord.out", output="f", plotfile="", 
          scale=1., fwhmpsf=2.5, interactive=no, radplots=no, verify=no, 
          update=no, emission=yes, sigma=INDEF, datamin=INDEF, 
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
          gcommands="", >& "dev$null")

          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
        mylist2 = "fff"

	while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, btakeout, npix )
}

#         print (btakeout, npix)
          btakeout = btakeout / npix

#	print (btakeout, radd, ix1, ix2, iy1, iy2)

#begin wht map correction 5/10/04

           delete ( "cord.out", verify-, >& "dev$null" )
           delete ( "f", verify-, >& "dev$null" )
           delete ( "ff.tab", verify-, >& "dev$null" )
           delete ( "fff", verify-, >& "dev$null" )

          print (xphot, yphot, > "cord.out")

           phot ("wcutout", "", 
           coords="cord.out", output="f", plotfile="", 
           scale=1., apertures=radd, calgorithm="none", zmag=25., 
           >& "dev$null")
          
          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
           mylist2 = "fff"
	   while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, wtakeout, npix )
}

          wtakeout = wtakeout / npix

#End the wht correction

# subtract BKGND mean from cutout -> ssi

          imarith ( cutout, "-", btakeout, ssi, title="", divzero=0., hparams="", pixtype="real", calctype="real", verbose-, noact- )


# compute initial rotation center

	  bgx = (x2 - x1) / 2 
	  bgy = (y2 - y1) / 2 

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

          radd = min( (ix2-ix1)/4, (iy2-iy1)/4 )
           
#     print (xcen, ycen, radd)
	if (xcen-radd-1 <= 0 || xcen+radd+1 > 4*radd) {cflag = yes 
                        flag = 1}
        if (ycen-radd-1 <= 0 || ycen+radd+1 > 4*radd) {cflag = yes
                        flag = 1}

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

           delete ( "cord.out", verify-, >& "dev$null" )
           delete ( "f", verify-, >& "dev$null" )
           delete ( "ff.tab", verify-, >& "dev$null" )
           delete ( "fff", verify-, >& "dev$null" )

          print (xphot, yphot, > "cord.out")

          phot (sqr, "", coords="cord.out", output="f", plotfile="", 
          scale=1., apertures=radd, calgorithm="none", zmag=25., >& "dev$null")

          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
           mylist2 = "fff"
	   while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, takeout, npix )
}
#problem
#	print (npix, radd)
           takeout = takeout / npix

#	print (takeout, npix)

# Takeout is the flux from the rotated and subtracted image	  

# Do the minimization: update best asymm value on each pass.

  	      if ( kk == 1 ) {
	        symm = takeout
	        centr = symm

   	        delete ( tmpdat, verify-, >& "dev$null" )
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
#	        
	   otakeout = takeout

#START

	if (x < 5) {x = x+50; flag = 1}
	if (y < 5) {y = y+50; flag = 1}

	  imget ( image, "i_naxis1" )
	  xdim = int(imget.value)
	  imget ( image, "i_naxis2" )
	  ydim = int(imget.value)

	  flprc; flprc

# get image dims

          size = xdim

       imdel ("fake.fits", yes, verify=no, default_acti=yes)

##8/4/03 - moved the mknoise of fake.fits inside the loop.  This is
## necessary such that it can be created with the correct dimentions.
## Below is the check that the edge is not being hit.

	  ix1 = int ( x - 49 )
          ix2 = int ( x + 50 )
	  iy1 = int ( y - 59 )
	  iy2 = int ( y + 50 )

#	print (ix1, ix2, iy1, iy2, xdim, ydim)

          xseg = 51
          yseg = 51

          if (ix1 <= 0) {ix1 = 1
                         xseg = x}
          if (ix2 <= 0) {ix2 = 3}
          if (iy1 <= 0) {iy1 = 1
                         yseg = y}
          if (iy2 <= 0) {iy2 = 3}


##8/4/03 - here is where the program checks to see if things are too
## high :
          if (ix1 >= xdim) {ix1 = xdim
                            xseg = xdim - x + 1}
          if (ix2 + dc >= xdim) {ix2 = xdim-1}
          if (iy1 + dc >= ydim) {iy1 = ydim-1}
          if (iy2 + dc >= ydim) {iy2 = ydim}

          xd = ix2 - ix1
          yd = iy2 - iy1

#here
#          print (btakeout)


          if (btakeout < 0) {btakeout = 0}


          mknoise ("fake", output="fake", title="", ncols=xd+1, nlines=yd+1,
          header="artdata$stdheader.dat", background=btakeout., gain=1., 
          rdnoise=stdd, poisson=yes, seed=1, cosrays="", ncosrays=0, 
          energy=30000., radius=0.5, ar=1., pa=0., comments=yes)

	xo = x
	yo = y

	flprc; flprc

	  section = "["//ix1//":"//ix2//","//iy1//":"//iy2//"]"
          segsec = "["//xseg//":"//xseg//","//yseg//":"//yseg//"]"

#   Make a (small) copy of the image to work with, called im
#	 	imdel ("wseg.fits", yes, verify=no, default_acti=yes)
	        imdel ("segZ.fits", yes, verify=no, default_acti=yes)
	        imdel ("segZ2.fits", yes, verify=no, default_acti=yes)
	        imdel ("segZ3.fits", yes, verify=no, default_acti=yes)
	     	imdel ("segZbg.fits", yes, verify=no, default_acti=yes)
	        imdel ("segZbg2.fits", yes, verify=no, default_acti=yes)
             	imdel ("segZbg3.fits", yes, verify=no, default_acti=yes)
	        imdel ("imbg.fits", yes, verify=no, default_acti=yes)
             
	        imdel ("ttest.fits", yes, verify=no, default_acti=yes)
	        imdel ("im.fits", yes, verify=no, default_acti=yes)
	        imdel ("im_pre.fits", yes, verify=no, default_acti=yes)

#print (section)
#         	imcopy ("test"//section, "ttest",  verbose- )
          	imcopy (image//section, "ttest",  verbose- )



#	print (section)
#SEGMENTATION MAP

# Shardha: 
#   Put a new segmentation map by doing  symbolic outside 
#   ln -s  /data/oslo3/jogee/gems/segmap/s1z02hot.seg.fits seg.fits
#
###SJ/JDbegin 
#	print (section)


        imcopy (segmap//section, "segZ", verbose-) 
###SJ/JDend	

        imcopy (wimage//section, "wseg", verbose-) 

   	  delete ( tmpdat, verify-, >& "dev$null" )

#Remove below if you don't know the num. of the galaxy, or name is not the SEG ID of the galaxy.

#print (segsec)
          imstat ( "segZ"//segsec, fields="max", lower=INDEF, upper=INDEF,  binwidth=0.1, format-, >> tmpdat )
	  mylist2 = tmpdat
	  jj = fscan ( mylist2, name) 

#         Store everything with 0 to -10  # Comment out to do pure new sky
          imreplace ("segZ", 427670., imaginary=0., lower=0., upper=0, radius=0.) # Comment this
 
#         Move all masks from 0 to name and name to max to 0
          imreplace ("segZ", 0., imaginary=0., lower=0., upper=name-1, radius=0.)
          imreplace ("segZ", 0., imaginary=0., lower=name+1, upper=327660, radius=0.) 

          imcopy ("segZ", "segZbg",  verbose- ) #for background computation, all objects at 1
 
#         Replace name and original 0 to 1, segZ2 is image with all other seg maps=0 
# Comment out 2nd one.
          imreplace ("segZ", 1., imaginary=0., lower=name, upper=name, radius=0.)
          imreplace ("segZ", 1., imaginary=0., lower=427670, upper=427670, radius=0.) # Comment this

          imreplace ("segZbg", 1., imaginary=0., lower=427670, upper=427670, radius=0.) 
          imreplace ("segZbg", 0., imaginary=0., lower=name, upper=name, radius=0.)

	  imarith ( "segZ", "*", "ttest", "segZ2", title="",  divzero=0.,
        	 hparams="", pixtype="", calctype="",verbose-, noact- )

	  imarith ( "segZbg", "*", "ttest", "segZbg2", title="",  divzero=0.,
        	 hparams="", pixtype="", calctype="",verbose-, noact- )

#         Set to 2 the galaxy/sky then 0 and 1 the other galaxies to 1
          imreplace ("segZ", 2., imaginary=0., lower=1, upper=2, radius=0.)
          imreplace ("segZ", 1., imaginary=0., lower=-1., upper=0, radius=0.)
          imreplace ("segZ", 0., imaginary=0., lower=2, upper=2, radius=0.)
	  
          imreplace ("segZbg", 2., imaginary=0., lower=1, upper=2, radius=0.)
          imreplace ("segZbg", 1., imaginary=0., lower=-1., upper=0, radius=0.)
          imreplace ("segZbg", 0., imaginary=0., lower=2, upper=2, radius=0.)
	  
#         Multiply this by the fake sky background and add together 
	  imarith ( "fake", "*", "segZ", "segZ3", title="",  divzero=0.,
        	 hparams="", pixtype="", calctype="",verbose-, noact- )
	  imarith ( "fake", "*", "segZbg", "segZbg3", title="",  divzero=0.,
        	 hparams="", pixtype="", calctype="",verbose-, noact- )

	  imarith ( "segZ3", "+", "segZ2", "im_pre", title="",  divzero=0.,
        	 hparams="", pixtype="", calctype="",verbose-, noact- )
	  imarith ( "segZbg3", "+", "segZbg2", "imbg", title="",  divzero=0.,
        	 hparams="", pixtype="", calctype="",verbose-, noact- )

           imstat ( "imbg", fields="mean", lower=INDEF, upper=INDEF,
        	   binwidth=0.1, format-, >> 'tmpdat' )

	  mylist2 = 'tmpdat'

	  jj = fscan ( mylist2, btakeout )
	  delete ( 'tmpdat', verify-, >& "dev$null" )


          imarith ( "im_pre", "-", btakeout, "im", title="",  divzero=0.,
        	    hparams="", pixtype="", calctype="",verbose-, noact- )


#finish computing im.

           delete ( "cntr", verify-, >& "dev$null" )
        imcntr ("im", xseg, yseg, cboxsize=15, > "cntr")

        mylist2 = "cntr"

	while ( fscan (mylist2, junk1, junk2, x, junk3, y) !=EOF )  
           delete ( "cntr", verify-, >& "dev$null" )

#again

        imcntr ("im", x, y, cboxsize=25, > "cntr")

        mylist2 = "cntr"

	while ( fscan (mylist2, junk1, junk2, x, junk3, y) !=EOF )  
           delete ( "cntr", verify-, >& "dev$null" )

#again

        imcntr ("im", x, y, cboxsize=5, > "cntr")

        mylist2 = "cntr"

	while ( fscan (mylist2, junk1, junk2, x, junk3, y) !=EOF )  
           delete ( "cntr", verify-, >& "dev$null" )

	xc = x
	yc = y

	kk = 2

	cflag = no
	iflag = no
	oflag = no
        halfcon = no

	tsb = 0

	  imget ( "im", "i_naxis1" )
	  xdim = int(imget.value)
	  imget ( "im", "i_naxis2" )
	  ydim = int(imget.value)

         size = min(xdim, ydim)

	  flprc; flprc
#print (x,y)

#   FIND THE CENTRAL SURFACE BRIGHTNESS

# 	BEGIN THE NEW SEARCH FOR THE RADII WHERE THE SURFACE BRIGHTENSS IS
# 	ONE HALF OF THE MAX.

	while ( !cflag ) {

	ix1 = int ( x - kk )
	ix2 = int ( x + kk )
	iy1 = int ( y - kk )
	iy2 = int ( y + kk )

        if ( x+kk+1 > xdim ) {cflag =yes
			      kk = xdim - x - 1
			      flag = 1}

        if ( y+kk+1 > ydim ) {cflag =yes
			      kk = ydim - y -1
			      flag = 1}

        if ( x-kk-1 <= 0 ) {cflag =yes
			   kk = x - 1
			   flag = 1}

        if ( y-kk-1 <= 0 ) {cflag =yes
	  		   kk = y - 1
			   flag = 1}
       
        if (kk > 198) {cflag = yes}

        section = "["//ix1//":"//ix2//","//iy1//":"//iy2//"]"

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
          scale=1., apertures=radd, calgorithm="none", zmag=25., >& "dev$null")

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
          scale=1., apertures=radd, calgorithm="none", zmag=25., >& "dev$null")

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

         if (osb < 1) {osb = 1; flag = 1; cflag = yes}
         if (npix < 1) {npix = 1; osb = 1; flag = 1; cflag = yes}

#	print (npix, nipix, osb, npix)
	tsb = ((osb - isb) / (npix - nipix)) / (osb/npix)

#       Remove lower # to print the Petrosian radii (tsb)
#	print ( kk, tsb, osb, tsb * npix / osb  )

#       Find the total radius, by the Petrosian radii at eta = 0.2

        t = tsb / (osb/npix)

#        print (tsb, radd)

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
 
#	print (btakeout, kk)

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

        if ( x+kk+1 > xdim ) {cflag =yes
			  kk = xdim - x - 1
			  flag = 1}

        if ( y+kk+1 > ydim ) {cflag =yes
			   kk = ydim - y -1
			   flag = 1}

        if ( x-kk-1 <= 0 ) {cflag =yes
			  kk = x - 1
			  flag = 1}

        if ( y-kk-1 <= 0 ) {cflag =yes
			   kk = y - 1
			   flag = 1}

	section = "["//ix1//":"//ix2//","//iy1//":"//iy2//"]"

           delete ( "cord.out", verify-, >& "dev$null" )
           delete ( "f", verify-, >& "dev$null" )
           delete ( "ff.tab", verify-, >& "dev$null" )
           delete ( "fff", verify-, >& "dev$null" )

          print (xphot, yphot, > "cord.out")
          radd = kk

          phot ("im", "", coords="cord.out", output="f", plotfile="", 
          scale=1., apertures=radd, calgorithm="none", zmag=25., >& "dev$null")

          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
         mylist2 = "fff"
	 while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, flux, npix )
}

#uncover the weight for the object

           delete ( "cord.out", verify-, >& "dev$null" )
           delete ( "f", verify-, >& "dev$null" )
           delete ( "ff.tab", verify-, >& "dev$null" )
           delete ( "fff", verify-, >& "dev$null" )

          print (xphot, yphot, > "cord.out")

          phot ("wseg", "", 
           coords="cord.out", output="f", plotfile="", 
           scale=1., apertures=radd, calgorithm="none", zmag=25., 
           >& "dev$null")
          
          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
           mylist2 = "fff"
	   while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, wobj, npix )
}
  

         if (npix < 1) {npix = 1; flag = 3; cflag = yes; iflag = yes; irad = 1; orad = 1}

   
           wobj = wobj / npix

#	print (wtakeout, wobj)

        if (wobj == 0) {wtakeout = 1}

	if (wtakeout > wobj) {wobj = wtakeout}


	takeout = takeout * sqrt(wtakeout/wobj)

#flux is the total flux of the galaxy within the eta = 0.2 radius

# 	print ('Begin computing concentrations')

	ll = 2

	tflux = 0

	while ( !iflag ) {

	ix1 = int ( x - ll )
	ix2 = int ( x + ll )
	iy1 = int ( y - ll )
	iy2 = int ( y + ll )

        if ( x+ll+1 > xdim ) {cflag =yes
			  ll = xdim - x - 1
			  flag = 1}

        if ( y+ll+1 > ydim ) {cflag =yes
			   ll = ydim - y -1
			   flag = 1}

        if ( x-ll-1 <= 0 ) {cflag =yes
			   ll = x - 1
			   flag = 1}

        if ( y-ll-1 <= 0 ) {cflag =yes
			   ll = y - 1
			   flag = 1}
 
        if (ll > 198) {cflag = yes
                       irad = 198
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
          scale=1., apertures=radd, calgorithm="none", zmag=25., >& "dev$null")

          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
         mylist2 = "fff" 
	 while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, tflux, npix )
}

#       Determines the inner radii for concentration, default = 0.2

        if ( tflux < -1) {
        tflux = sqrt(tflux*tflux)
        flag = 1}


	if (tflux < 0 || tflux == 0) {
        iflag = yes
	irad = 1
        flag = flag + 1}

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

        if ( x+ll+1 > xdim ) {cflag =yes
			    ll = xdim - x - 1
			    flag = 1}

        if ( y+ll+1 > ydim ) {cflag =yes
			    ll = ydim - y -1
			    flag = 1}

        if ( x-ll-1 <= 0 ) {cflag =yes
			   ll = x - 1
			   flag = 1}

        if ( y-ll-1 <= 0 ) {cflag =yes
			   ll = y - 1
			   flag = 1}

        if (ll > 198) {cflag = yes
                       oflag = yes
                       orad = 198}

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
         scale=1., apertures=radd, calgorithm="none", zmag=25., >& "dev$null")

          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
         mylist2 = "fff"
	 while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, tflux, npix )
}

#       Determine half-light radius


        if (tflux <= 0) {tflux = 10
                         otflux = tflux + 1
                         flag = flag + 1
                         oflag = yes}


	if (tflux < 0 || tflux == 0) {
        halfcon = yes
        flag = flag + 1}

        if (!halfcon) {

        if ( tflux > 0.5 * flux) {
          halfcon = yes

          halfr = ll - ((0.5*flux - tflux)/(otflux - tflux))}}

# Deteremines the outer radii for concentration, default = 0.8

	if (tflux < 0 || tflux == 0) {
        oflag = yes
        flag = flag + 1
	tflux = flux}
!
	if ( tflux > 0.8 * flux) {
	oflag = yes
	orad = ll

        orad = ll - ((0.8*flux - tflux)/(otflux - tflux))}

	if ( ll == 0 ) {
	oflag = yes}

	ll = ll + 1 }

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

#	print ('Begin Computation of A')

#############################################################################
#                                                                           #
#
#                      Begin Asymmetry Source computation
#                        Note 'takeout' computed earlier                    #
#                                                                           #
#############################################################################
       	  imget ( "im", "i_naxis1" )
	  xdim = int(imget.value)
	  imget ( "im", "i_naxis2" )
	  ydim = int(imget.value)

	  flprc; flprc

# get image dims

          size = min(xdim, ydim)
  
          rad = bigrad / 1.5
	  rad0 = rad

#          print (x, rad0, rad, dc, size, ydim)
	  if ( x+rad+dc+1>size || x-rad-dc-1<0 || y+rad+dc+1>size || y-rad-dc-1<0 ) {
	     if ( x+rad0+dc+1>size ) rad = min ( rad-dc-1, xdim-x-dc-1 )
	     if ( y+rad0+dc+1>size ) rad = min ( rad-dc-1, ydim-y-dc-1 )
	     if ( x-rad0-dc-1<0 )    rad = min ( rad-dc-1, x-dc-1 )
	     if ( y-rad0-dc-1<0 )    rad = min ( rad-dc-1, y-dc-1 ) 
	if (rad < 0) rad = 0
	     print ( 'WARNING: user radius is off image; ' )
	     print ( 'radius reset from ',rad0,' to ',rad )
	
	  }
#          print (x, rad0, rad, dc, size, ydim)

           if (rad < 1) {rad = 1
                         dc = 1}


# define extraction box and extract

	  xorg = x
	  yorg = y

	  ixx1 = int ( xorg-rad)
          ixx2 = int ( xorg+rad)

          iyy1 = int ( yorg-rad)
          iyy2 = int ( yorg+rad)

#         print (x,y,xorg,yorg,rad)

	  psect = "["//ixx1-dc//":"//ixx2+dc//","//iyy1-dc//":"//iyy2+dc//"]"
	  section = "["//dc//":"//ixx2-ixx1+dc//","//dc//":"//iyy2-iyy1+dc//"]"
#	print (psect)
          imcopy ("im"//psect, ssi,  verbose- )

# subtract BKGND

#          imarith ( cutout, "-", btakeout, ssi, title="",  divzero=0.,
#        	    hparams="", pixtype="", calctype="",verbose-, noact- )

# create abs value image and get statistics for denominator value

	  imfunction ( ssi, imsqr, "abs", verbose- )
#          imcopy (ssi, imsqr,  verbose- )

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
			
        if ( xphot-radd-1 <= 0 ) {cflag =yes
			   radd = rad + dc - 1
			   flag = flag + 1}

        if ( yphot-radd-1 <= 0 ) {cflag =yes
			   radd = rad + dc - 1
			   flag = flag + 1}

          phot (imsqr, "", coords="cord.out", output="f", plotfile="", 
          scale=1., apertures=radd, calgorithm="none", zmag=25., >& "dev$null")
          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
         mylist2 = "fff"
	 while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, bottom, npix )
}

         if (npix < 1) {npix = 1; flag = 1; cflag = yes}


           bottom = bottom / npix

  	        delete ( tmpdat, verify-, >& "dev$null" )

          imstat ( imsqr//section, fields="stddev", lower=INDEF, upper=INDEF,
        	   binwidth=0.1, format-, >> tmpdat )

	  mylist2 = tmpdat

	  jj = fscan ( mylist2, istd )
	  delete ( tmpdat, verify-, >& "dev$null" )
	  imdelete ( imsqr, verify-, >& "dev$null" )

# compute initial rotation center

	  x = rad + dc
	  y = rad + dc

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
		       yout=ycen, ncols=INDEF, nlines=INDEF, interpolant="linear", boundary="nearest", constant=0., nxblock=256, nyblock=256, 
verbose- )

              imarith ( rot, "-", ssi, sub, title="", divzero=0.,
		        hparams="", pixtype="real", calctype="real", verbose-, noact- )
 	      imfunction ( sub, sqr, "abs", verbose- )

# get mean value of abs. val. of resid. image: top

           delete ( "cord.out", verify-, >& "dev$null" )
           delete ( "f", verify-, >& "dev$null" )
           delete ( "ff.tab", verify-, >& "dev$null" )
           delete ( "fff", verify-, >& "dev$null" )
	   
#	   print (xcen, ycen, radd)

        if ( xcen+radd+1 > 2*rad + 2*dc ) {cflag =yes
			    radd = 2*rad + 2*dc - xcen
			    flag = flag + 1}

        if ( ycen+radd+1 > 2*rad + 2*dc ) {cflag =yes
			    radd = 2*rad + 2*dc - ycen
			    flag = flag + 1}
			
        if ( xcen-radd-1 <= 0 ) {cflag =yes
			   radd = xcen - 1
			   flag = flag + 1}

        if ( ycen-radd-1 <= 0 ) {cflag =yes
			   radd = ycen - 1
			   flag = flag + 1}
          xphot = xcen
 	  yphot = ycen

          print (xphot, yphot, > "cord.out")

          phot (sqr, "", coords="cord.out", output="f", plotfile="", 
          scale=1., apertures=radd, calgorithm="none", zmag=25., >& "dev$null")
          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
         mylist2 = "fff"
	 while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, top, npix )
}

#        print (top, npix, radd)

        if (npix < 1) {npix = 1; flag = 1; cflag = yes}

         top = top / npix

# Do the minimization: update best asymm value on each pass.
 
      if (bottom < 1) {bottom = 1; flag = 5; cflag = yes}

	      if ( kk == 1 ) {
	        symm = (top - takeout) / (bottom)
	        centr = symm
	      }

	      tsymm = (top - takeout) / (bottom)

#	print (tsymm, radd)

	      if ( tsymm < symm ) {
#	print (top, takeout, bottom)
	        symm = tsymm
	        x = xcen
	        y = ycen}
#	print (top, takeout, bottom)           	      

# Error computation

	      error = (stdd / ((bottom+istd)))

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

#	  print ('Begin Computation of S')

	flprc
	flprc


#	print ('begin symm')

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

          imcopy ("im", "temp1",  verbose- )

	  boxcar ("temp1", ssi, 5, 5, boundary="nearest", constant=0.)

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

        if ( x-radd-1 <= 0 ) {cflag =yes
			      radd = x - 1
			      flag = flag + 1}

        if ( y-radd-1 <= 0 ) {cflag =yes
			      radd = y -1
			      flag = flag + 1}

          phot ("c1", "", coords="cord.out", output="f", plotfile="", 
          scale=1., apertures=radd, calgorithm="none", zmag=25., >& "dev$null")

          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
         mylist2 = "fff"
	 while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, top, npix )
}

	if (npix < 1) {npix =1; cflag = yes; flag = 3}

           top = top / npix

#Flux from the total image

           delete ( "cord.out", verify-, >& "dev$null" )
           delete ( "f", verify-, >& "dev$null" )
           delete ( "ff.tab", verify-, >& "dev$null" )
           delete ( "fff", verify-, >& "dev$null" )

          print (xphot, yphot, > "cord.out")

          phot (ssi, "", coords="cord.out", output="f", plotfile="", 
          scale=1., apertures=radd, calgorithm="none", zmag=25., >& "dev$null")

          tae ("f", "ff", 1, "27", "", datatype="", colunits="", colfmt="")
          tae ("f", "ff", 1, "28", "", datatype="", colunits="", colfmt="")
     
          tprint ("ff", prparam=no, prdata=yes, pwidth=80, plength=0, 
          showrow=yes, orig_row=yes, showhdr=yes, showunits=yes, columns="", 
          rows="-", option="plain", align=yes, sp_col="", lgroup=0, > "fff")
      
         mylist2 = "fff"
	 while ( fscan (mylist2, junk1) !=EOF )   {
	   ii = fscan ( mylist2, junk2, bottom, npix )
}

	if (npix < 1) {npix = 1; cflag = yes; flag = 3}

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

          if (ix1 - dc <= 0) {ix1 = 1 + dc}
          if (ix2 - dc <= 0) {ix2 = 3 + dc}
          if (iy1 - dc <= 0) {iy1 = 1 + dc}
          if (iy2 - dc <= 0) {iy2 = 3 + dc}

          if (ix1 + dc >= 0) {ix1 = xdim - 3 - dc}
          if (ix2 + dc >= 0) {ix2 = xdim - 1 - dc}
          if (iy1 + dc >= 0) {iy1 = ydim - 3 - dc}
          if (iy2 + dc >= 0) {iy2 = ydim - 3 - dc}

		bsect = "["//ix1-dc//":"//ix2+dc//","//iy1-dc//":"//iy2+dc//"]"
          
#print (bsect)

        	imcopy (image//bsect, "cutout",  verbose- )

          	imarith ("cutout", "-", btakeout, "temp2", title="",  divzero=0., hparams="", pixtype="", calctype="",verbose-, noact- )

	  	boxcar ("temp2", ssi, 5, 5, boundary="nearest", constant=0.)

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

        if ( xphot-radd-1 <= 0 ) {cflag =yes
			      radd = xphot - 1
			      }

        if ( yphot-radd-1 <= 0 ) {cflag =yes
			      radd = yphot - 1
			      }

          phot ("c1", "", coords="cord.out", output="f", plotfile="", 
          scale=1., apertures=radd, calgorithm="none", zmag=25., >& "dev$null")

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
    	   ttakeout = ttakeout * sqrt(wtakeout/wobj)

#Delete files

          imdelete ( "cutout", verify-, >& "dev$null" )
          imdelete ( cutout, verify-, >& "dev$null" )
          imdelete ( ssi, verify-, >& "dev$null" )
          imdelete ( rot, verify-, >& "dev$null" )
          imdelete ( "c1", verify-, >& "dev$null" )
          imdelete ( "t1", verify-, >& "dev$null" )
          imdelete ( "c2", verify-, >& "dev$null" )
          imdelete ( "temp1", verify-, >& "dev$null" )
          imdelete ( "temp2", verify-, >& "dev$null" )

# Reformat s and error to take 3 sig. dig., and output

#	print (top, ttakeout, bottom)

	if (bottom < 0) {bottom = 1; flag = 3; cflag = yes}

          s = (top - ttakeout) / bottom 
	  se = 10*(sqrt(abs(npix*(top-ttakeout))))/(npix*(bottom+istd))
	  s = 10 * (real(int(1000. * s))/1000.)	  
	  se = real(int(1000. * se))/1000. 

flprc
flprc
flprc

#################################################################
#Start the M_20/Gini stuff
#################################################################

	radd = radd * 1.5


#print (x, y, rad, dc, size)

        if ( x+rad+dc > size ) {
			      rad = size - x - dc
			      flag = flag + 1}

        if ( y+rad+dc > size ) {
			      rad = size - y -dc
			      flag = flag + 1}

        if ( x-rad-dc < 1 ) {
			      rad = x - 1 - dc
			      flag = flag + 1}

        if ( y-rad-dc < 1 ) {
			      rad = y - 1 - dc
			      flag = flag + 1}
        if (( x+rad+dc < 1 ) || ( y+rad+dc < 1 ) || ( x-rad-dc < 1 ) || ( y-rad-dc < 1 )) {rad = 0; x = 5; dc = 2; flag = 999}

#print (x, rad, dc, size)
#stop


	  section = "["//int(x-rad-dc)//":"//int(x+rad+dc)//","//int(y-rad-dc)//":"//int(y+rad+dc)//"]"

#print (section, rad)

	listpixels ("im"//section, wcs="logical", formats="", verbose=no, >> 'tab')

        mylist2 = "tab2"

	xm = 1+rad+dc
	ym = 1+rad+dc

	sort ("tab", column=3, ignore_white=no, numeric_sort=yes, reverse_sort=yes, >> 'tab2')

	fdis = 0
	tdis = 0
	tabflux = 0
#GINI
       numc = 1

	while (fscan(mylist2, xt, yt, inten) !=EOF)
        {

	dis = ((xm-xt)**2 + (ym-yt)**2)
	
 	if (sqrt(dis) < rad) {fdis = fdis + inten*dis
        numc = numc+1
	tabflux = tabflux + inten }

	   if (tabflux < 0.2*flux && sqrt(dis) < rad) {tdis = tdis + inten*dis}
	print (x, xt, y, yt, rad, tabflux, flux, fdis, tdis, inten, dis, >> 'test.out')

 }             

#GINI  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
       aveflux = tabflux / numc
       i = numc+1
       gini2 = 0

       if (numc == 1) numc = 2

       mylist2 = "tab2"

       while (fscan(mylist2, xt, yt, inten) !=EOF)
       {       dis = ((xm-xt)**2 + (ym-yt)**2)

       if (sqrt(dis) < rad)  {i = i - 1
         gini1 = ((2*i) - numc - 1)*inten
         gini2 = gini2 + gini1
#          print (i, gini1, gini2)
}}

        if (aveflux == 0) aveflux = 1

         gini = (1/(aveflux*numc*(numc-1))) * gini2
#GINI >>>>>>>>>>>>>>>>>>>

#          print (aveflux, numc, gini2)

           delete ( "tab", verify-, >& "dev$null" )
           delete ( "tab2", verify-, >& "dev$null" )

#	print (tdis, fdis)

	if (tdis <= 0 || fdis <= 0) {tdis = fdis*fdis+1; fdis = tdis}

         if (fdis < 1) set fdis = 1

#	print (tdis, fdis)

	c20 = log10(tdis/fdis)

	  c20 = real(int(1000. * c20))/1000. 
         gini = real(int(1000.*gini))/1000

flprc
flprc
flprc



#####################
#END OF C_20
#####################

#END the computation
 
	delete ( tmpdat, verify-, >& "dev$null" )

# Reformat symm and error to take 3 sig. dig., and output

	  symm = real(int(1000. * symm))/1000.
	  error = real(int(1000. * error/3))/1000
	  xorg = real(int(1000. * xorg))/1000
	  yorg = real(int(1000. * yorg))/1000
	  xp = real(int(1000. * (xp+xo-xc)))/1000
          yp = real(int(1000. * (yp+yo-yc)))/1000 
	  rad = real(int(1000. * rad))/1000
	
	  conc = real(int(1000. * conc))/1000.
	  econc = real(int(1000. * econc/3))/1000
	  s = real(int(1000. * s))/1000
	  se = real(int(1000. * se))/1000
	  bigrad = real(int(1000. * bigrad))/1000

#renormalize the flux

	  flux = flux / 12000

	  flux = real(int(1000. * flux))/1000
	  be = bigrad
	  be = real(int(1000. * be))/1000.

   	printf ("%10s %8f %8f %5f %5f %5f %5f %5f %5f %5f %6f %6f %6f %6f %6f %6f\n", id, xp, yp, conc, econc, symm, error, s, se, c20, gini, bigrad, 
halfr, flux, sflag, flag, >> ofile)

	print  (id, " ", conc, econc, symm, error, s, se, c20, gini, bigrad, halfr, flux, sflag, flag )

	flag = 0

        flprc
        flprc

        dc = 5

	takeout = otakeout
 
 }
end
 
#display ("im",
#1, bpmask="", bpdisplay="none", bpcolors="red", overlay="",
#ocolors="green", erase=yes, border_erase=no, select_frame=yes, repeat=no,
#fill=no, zscale=yes, contrast=0.25, zrange=yes, zmask="", nsample=1000,
#xcenter=0.5, ycenter=0.5, xsize=1., ysize=1., xmag=1., ymag=1., order=0,
#z1=0., z2=1., ztrans="log", lutfile="")

#imexam ("im",
#1, "", output="", ncoutput=101, nloutput=101, logfile="", keeplog=no,
#defkey="a", autoredraw=yes, allframes=yes, nframes=0, ncstat=5, nlstat=5,
#graphcur="", imagecur="", wcs="logical", xformat="", yformat="",
#graphics="stdgraph", display="display(image='$1',frame=$2)", use_display=yes)

 

