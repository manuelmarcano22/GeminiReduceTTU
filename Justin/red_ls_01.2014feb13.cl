# GMOS LongSlit reduction script - part 1 of 3
# This part does the wavelength independent pre-reduction (bias, dark, arcs calibration) 
# ------------------------------------------------------------------------------------
# GG, Jan 2014

# ** Assumes all necessary files (bias,flats,arcs and science) are in a common, local directory **
# I suggest to use a separate directory per target for cleanliness
# Builds relevant lists using hselect

gemini
gmos
delete.verify=no

###############################################################
print "----------------Building lists-----------------"

string *list1

delete Flat.lst,Flat.txt,Dark.lst,Dark.txt,Bias.lst,Bias.txt,Sci.lst,Sci.txt,Arc.lst,Arc.txt,wavelengths.txt verify-
delete gprepare_arc.log,gsreduce_arc.log,gswav_arc.log verify-

delete gprepare_sci.log
delete gsreduce_sci.log
unlearn gprepare
unlearn gsreduce

# add the mdf otherwise we'll have trouble with gscut later on!!
gprepare.fl_addmdf=yes
gprepare.logfile="gprepare_sci.log"

hselect S*.fits[0] $I "OBJECT='CuAr'" > Arc.lst
! awk '{gsub(/.fits\[0\]/,""); print}' Arc.lst > Arc.txt

hselect S*.fits[0] $I "OBJECT='GCALflat'" > Flat.lst
! awk '{gsub(/.fits\[0\]/,""); print}' Flat.lst > Flat.txt

hselect S*.fits[0] $I "OBJECT='Bias'" > Bias.lst
! awk '{gsub(/.fits\[0\]/,""); print}' Bias.lst > Bias.txt

## Uncomment these if using for science
#hselect S*.fits[0] $I "OBSCLASS='science' && OBSTYPE='OBJECT'" > Sci.lst
#! awk '{gsub(/.fits\[0\]/,""); print}' Sci.lst > Sci.txt

## Uncomment these if using for standard stars 
hselect S*.fits[0] $I "OBSCLASS='partnerCal' && OBSTYPE='OBJECT'" > Sci.lst
! awk '{gsub(/.fits\[0\]/,""); print}' Sci.lst > Sci.txt

hselect @Arc.lst waveleng yes >> wavelengths.txt

###############################################################
print "----------------Bias reduction-----------------"
gbias (inimages="@Bias.txt", outbias="bias_out", fl_over+, fl_trim+, fl_vardq-, rawpath="./", logfile=gbias.log)

#imstat bias_out.fits[2][200:600,2000:2400]




###############################################################
print "----------------GCal Flats reduction-----------------"
# non interactive works well, order 7

list1='Flat.txt'
while (fscan (list1, s1) != EOF) {
	gprepare.fl_addmdf=yes
	print(s1)
# DOESN'T WORK IF YOU USE ='yes' instead of '+'
	gsflat (s1,"f"//s1, bias="bias_out", fl_fixpix+, fl_inter-, function='spline3', order='19', fl_detec+, fl_double-, nshuffle=1536, logfile="gsflat.log")
}
;
ls fS*.fits >> Flat_red.txt

###############################################################
print "----------------Arcs reduction-----------------"

gprepare @Arc.txt outpref="p" rawpath="./" fl_addmdf+ 

delete Arc.txt
ls pS*.fits >> Arc.txt

gsreduce @Arc.txt outpref="u" logfile="gsreduce_arc.log" rawpath="./" fl_over+ fl_trim+ fl_bias+ bias="bias_out" fl_dark- fl_flat- fl_gmosaic+ fl_fixpix+ fl_gsappwave+ fl_cut+ ovs_flinter- fl_vardq-

# Wavelenght calibration
# order 6 works OK (0.1 Angstrom error). Order 7 not so well. higher orders - higher nonlinear term. 
gswavelength up* fl_inter- order=6 nlost=10 logfile="gswav_arc.log"

ls upS*.fits >> Arc_red.txt

list1 = 'Arc_red.txt'
while (fscan (list1, s1) != EOF) {
	gstransform(s1,wavtran=s1,logfile="gstransform_arc.log")
}




###############################################################
print "----------------Science pre-reduction-----------------"
# Science: first run of gsreduce for the wavelength-indepentant stage (bias and dark)
gsreduce @Sci.txt outpref="gsr" logfile=gsreduce_sci.log rawpath=./ fl_over+ fl_trim+ \
fl_bias+ bias=bias_out fl_dark- fl_flat- fl_gmosaic- fl_fixpix- fl_gsappwave- fl_cut- ovs_flinter- fl_vardq-

delete Sci_red.txt
ls gsr*.fits > Sci_red.txt

print "run 'gemextn' on the reduced science files and verify that the MDF extension is there"
print "Now you should have the lists Sci_red.txt, Flat_red.txt and Arc_red.txt"
print "Verify that they exist and are non-empty"
ls -ltr *.txt

