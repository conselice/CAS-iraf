procedure sort4 

struct *mylist1		{ prompt = "wiz" }
struct *mylist2
struct *mylist10
begin


        real id, ra, dec, x, y, Bmag, eBmag, Vmag, eVmag, Imag, eImag, Zmag, eZmag, Jmag, eJmag, Hmag, eHmag, ZB, ZBmin, ZBmax, TB, ODDs, Chisq
        real takeout, npix, new, dc, maxval
        real xc, yc
        int x1, x2, y1, y2, xx1, xx2, yy1, yy2
        int ii, xdim, ydim, out, num, jj
        int bgx1, bgx2, bgy1, bgy2
        real b, be, r, re, i, ie, pz, bk, rk, ke, ik, je
        real kim, field
        string junk2, junk1, name, namee, image, segi, whti, simage, sect, wht, seg, sect2


stsdas
imgtools

############################################################################
# Set up
############################################################################

# make and clear temp images and fi
# get inputs

	mylist1 = "list"
	mylist10 = "cata"


        ii = 0
        dc = 50

#changed dc from 150 to 200 then to 50 for wfc3

#changed	while ( fscan ( mylist10, new, ra, dec) != EOF) {

	while ( fscan ( mylist10, id, ra, dec) != EOF) {

        ii = ii + 1 

        num = 0
        xc = 51
        yc = 51

#
	mylist1 = "list"

	while ( fscan ( mylist1, image, seg, wht) != EOF) {
  
          imget ( image, "i_naxis1" )
          xdim = int(imget.value) 
          imget ( image, "i_naxis2" )
          ydim = int(imget.value)

           rd2xy (image, ra, dec, hour=no, >> "out")

	print (id, image, ra, dec)

          mylist2 = "out"
           while ( fscan (mylist2, junk1, junk2,x,junk2,junk1, y) !=EOF )
          delete ( "out", verify-, >& "dev$null" ) 


          if (x > 3 && x+3 < xdim && y > 3 && y+3 < ydim)  {
                 x = x 
                 y = y 
                 namee = image  

          num = num + 1

          x1 = int (x-dc)
          x2 = int (x+dc)
          y1 = int (y-dc)
          y2 = int (y+dc)

          xx1 = int (x-2)
          xx2 = int (x+2)
          yy1 = int (y-2)
          yy2 = int (y+2)

print (xdim, ydim, x, y)

          if (x1 < 1) {xc = xc + x1
                        x1 = 1}
          if (x2 > xdim) {x2 = xdim-1}

          if (y1 < 1) {yc = yc + y1; 
                         y1 = 1}
          if (y2 > ydim) {y2 = ydim-1}

          out = int (id)

           sect = "["//x1//":"//x2//","//y1//":"//y2//"]"
           sect2 = "["//xx1//":"//xx2//","//yy2//":"//yy2//"]"



               imstat (image//sect2, fields="max", lower=INDEF,
		         upper=INDEF, binwidth=0.1, format-, >> 'tmpdat' )

                mylist2 = 'tmpdat'
                jj = fscan ( mylist2, maxval )
     	        delete ('tmpdat', verify-, >& "dev$null" )

	      if ( maxval > 0 ) {

 #changed below for all images created to have a u. in front
          imcopy (image//sect, "u."//out//"."//num, verbose=yes)
          imcopy (wht//sect, "u.wht."//out//"."//num, verbose=yes)
          imcopy (seg//sect, "u.seg."//out//"."//num, verbose=yes)

           print ("u."//out//"."//num, ' ', id, xc, yc, "u.wht."//out//"."//num, ' ', "u.seg."//out//"."//num, ' ', 0, >> 'input.cas')

	print (out//"."//num, ' ', id, ra, dec, >> 'catalog.cas') }}

	      		    
}}
         
flprc
flprc


end





#ncas ("egs.input", "output")

