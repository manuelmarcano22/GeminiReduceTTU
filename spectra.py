# For reduced:
from stsci.tools import capable
capable.OF_GRAPHICS = False

from pyraf import iraf
from pyraf.iraf import gemini, gemtools, gmos, onedspec
import fileSelect as fs
import copy, os
from astropy.io import fits


#For plotting
from astropy.io import fits
import os
from stsci.tools import capable
capable.OF_GRAPHICS = False
from pyraf import iraf
import numpy as np
from shutil import copyfile

import pyds9


#Bokeh plotting
#from bokeh.io import push_notebook, show, output_notebook
from bokeh.plotting import Figure
from bokeh.models import Span, Label, Arrow, NormalHead
from bokeh.models import HoverTool, tools, ColumnDataSource, CustomJS, Slider, BoxAnnotation
from bokeh.layouts import  column, row
from bokeh.palettes import viridis, BrBG8, Accent8, OrRd8, Pastel2_4, Category20b_8, Viridis8,Dark2_8,Paired8, Set1_4
import re

from bokeh.plotting import output_notebook, figure, show
import numpy as np

from math import floor



output_notebook()



name = pyds9.ds9_targets()
if name != 'None':
    try:
        name1 = name[0].split()[1]
        #print(name1)
        d = pyds9.DS9(name1)
    except:
        pass
    
    

##Emission Lines
class line(object):
    def __init__(self,name):
        self.name = name
        
def reducedgemini(qd,arcname,datadirall='../raw/',ds9=True,skyreg='100:500'):
    
    #path
    databasename = datadirall+'obsLog.sqlite3'
    dbFile= datadirall+'obsLog.sqlite3'

    
    gemtools.gemextn.unlearn()    # Disarm a bug in gbias
    gmos.gbias.unlearn()
    biasFlags = {
        'logfile':'biasLog.txt','rawpath':'../raw/','fl_vardq':'yes','verbose':'no'
    }
    regions = ['Full']
    for r in regions:
        # The following SQL generates the list of full-frame files to process.
        SQL = fs.createQuery('bias', qd[r])
        biasFiles = fs.fileListQuery(dbFile, SQL, qd[r])

        # The str.join() funciton is needed to transform a python list into a
        # comma-separated string of file names that IRAF can understand.
        if len(biasFiles) > 1:
            gmos.gbias(','.join(str(x) for x in biasFiles), 'MCbias'+r,
                       **biasFlags)
            
    #Clean up after BIAs

    iraf.imdel('gS2018*.fits')
    
    # Open it in DS9 to see it
    if ds9 == True:
        try:
            biasfit = fits.open('MCbiasFull.fits')
            d.set_fits(biasfit)
        except:
            pass
    
        # Set the task parameters.
    qd['Full'].update({'DateObs':'*'})
    gmos.gireduce.unlearn()
    gmos.gsflat.unlearn()
    # The response fitting should be done interactively.
    flatFlags = {
        'fl_over':'yes','fl_trim':'yes','fl_bias':'yes','fl_dark':'no',
        'fl_fixpix':'no','fl_oversize':'no','fl_vardq':'yes','fl_fulldq':'yes',
        'rawpath':'../raw/','fl_inter':'no','fl_detec':'yes',
        'function':'spline3','order':'13,11,28',
        'logfile':'gsflatLog.txt','verbose':'no'
        }
    for r in regions:
        qr = qd[r]
        flatFiles = fs.fileListQuery(dbFile, fs.createQuery('gcalFlat', qr), qr)
        if len(flatFiles) > 0:
            gmos.gsflat (','.join(str(x) for x in flatFiles), 'MCflat'+r,
                    bias='MCbias'+r, **flatFlags)

    iraf.imdel('gS2018*.fits,gsS2018*.fits')

    if ds9 == True:
        try:
            flatfile = fits.open('MCflatFull.fits')
            d.set_fits(flatfile)
        except:
            pass
        
        
    # Set task parameters.
    gmos.gsreduce.unlearn()
    sciFlags = {
        'fl_over':'yes','fl_trim':'yes','fl_bias':'yes','fl_gscrrej':'no',
        'fl_dark':'no','fl_flat':'yes','fl_gmosaic':'yes','fl_fixpix':'no',
        'fl_gsappwave':'yes','fl_oversize':'no',
        'fl_vardq':'yes','fl_fulldq':'yes','rawpath':'../raw/',
        'fl_inter':'no','logfile':'gsreduceLog.txt','verbose':'no'
    }
    arcFlags = copy.deepcopy(sciFlags)
    arcFlags.update({'fl_flat':'no','fl_vardq':'no','fl_fulldq':'no'})


    # Arc exposures
    for r in regions:
        qr = qd[r]
        arcFiles = fs.fileListQuery(dbFile, fs.createQuery('arc', qr), qr)
        if len(arcFiles) > 0:
            gmos.gsreduce (','.join(str(x) for x in arcFiles), bias='MCbias'+r,
                      **arcFlags)


    # Science exposures
    r = 'Full'
    sciFiles = fs.fileListQuery(dbFile, fs.createQuery('sciSpec', qd[r]), qd[r])
    if len(sciFiles) > 0:
        gmos.gsreduce (','.join(str(x) for x in sciFiles), bias='MCbias'+r,
                  flatim='MCflat'+r, **sciFlags)
        
    iraf.imdel('gS2018*.fits')
    
    
    #############Arccc
        # Set task parameters
    gmos.gswavelength.unlearn()
    waveFlags = {
        'coordlist':'gmos$data/CuAr_GMOS.dat','fwidth':6,'nsum':50,
        'function':'chebyshev','order':5,
        'fl_inter':'no','logfile':'gswaveLog.txt','verbose':'no'
        }
    # Must select specific wavecals to match science exposures.
    #prefix = 'gsS20180708S00'
    prefix = 'gs'

    #for arc in [arcname+'.fits']:
    for arc in [arcname]:
         #gmos.gswavelength (prefix+arc, **waveFlags)
        gmos.gswavelength (prefix+arc, **waveFlags)
        
        
        
    # Set task parameters.
    gemtools.gemcombine.unlearn()
    sciCombFlags = {
        'combine':'average','reject':'ccdclip',
        'fl_vardq':'yes','fl_dqprop':'yes',
        'logfile':'gemcombineLog.txt.txt','verbose':'no'
    }
    stdCombFlags = copy.deepcopy(sciCombFlags)
    stdCombFlags.update({'fl_vardq':'no','fl_dqprop':'no'})
    gmos.gstransform.unlearn()
    transFlags = {
        'fl_vardq':'yes','interptype':'linear','fl_flux':'yes',
        'logfile':'gstransLog.txt'
    }
    # The sky regions should be selected with care, using e.g. prows/pcols:
    #   pcols ("tAM2306b.fits[SCI]", 1100, 2040, wy1=40, wy2=320)
    gmos.gsskysub.unlearn()
    skyFlags = {
        'fl_oversize':'no','fl_vardq':'yes','logfile':'gsskysubLog.txt'
    }        
    
    
    
    sciTargets = {
    qd['Full']['Object']:{'arc':'gs'+arcname,'sky':skyreg},
    }

    prefix = "gs"

    for targ,p in sciTargets.iteritems():
        qs = qd['Full']
        qs['Object'] = targ
        # Fix up the target name for the output file
        sciOut = targ.split('-')[0]+targ[-1]
        sciFiles = fs.fileListQuery(dbFile, fs.createQuery('sciSpec', qs), qs)
        gemtools.gemcombine (','.join(prefix+str(x) for x in sciFiles),
                             sciOut, **sciCombFlags)
        gmos.gstransform(sciOut, wavtraname=p['arc'], **transFlags)
        gmos.gsskysub('t'+sciOut, long_sample=p['sky'], **skyFlags)

    iraf.imdel("gsS2018*.fits")

    if ds9 == True:
        try:
            spectrafile = fits.open('t'+sciOut+'.fits')
            d.set_fits(spectrafile)
        except:
            pass
        


    
def plotapertures(sourcenameorg):  
    
    #Databse dowlaoded from old agithub vimos cx25
    databasegeneral = 'database/apcx25sexm'
    filename = 'database/ap'+sourcenameorg
    sourcename = 'copy'+sourcenameorg
    
    if not os.path.isfile(sourcename+'.fits'):
        iraf.imcopy(sourcenameorg+'[SCI]',sourcename)
    
    apefile = fits.open(sourcename+'.fits')
    ##For srfm[0].header["CTYPE1"] = 'LINEAR'
    xn = apefile[0].header["NAXIS1"]
    refx = apefile[0].header["CRVAL1"]
    #step = apefile[2].header['CD1_1']
    step = int(apefile[0].header['CD1_1'])
    #Had to do this since it was given me an equalt stepo size. 
    cr = apefile[0].header['CRPIX1']
    ape = apefile[0].data
    #
    xlist = [ refx + step*(i - cr) for i in np.arange(ape.shape[1] +1) ]
    #

    #Default is half
    y = ape[:,int(ape.shape[1]/2)]
    y = ape[:,0]

    x = list(range(y.shape[0]))
    c2 = Figure(x_axis_label='index',title='Wavelength',toolbar_location="above",active_scroll='wheel_zoom')
    source2 = ColumnDataSource(data=dict(x=x,y=y))

    #Separation to plot. Plot eevery given wavelength
    sepaape = 10


    sourceall2 = ColumnDataSource(data= dict([ (str(int(xlist[i])), ape[:,i]  ) for i in list(range(0,ape.shape[1],sepaape)) ]))

    wavelist =  [ int(xlist[i]) for i in list(range(0,ape.shape[1],sepaape)) ]
    sepawave = wavelist[1]- wavelist[0]
    #sepawave = sepaape

    c2.line('x','y', source=source2)

    callback2 = CustomJS(args=dict(source2 = source2, sourceall2 = sourceall2 ), code="""

            var data2 = source2.data;
            var data22 = sourceall2.data;
            var f = cb_obj.value;
            y = data2['y'];
            y2 = data22[f.toString()];


            for (i = 0; i < y.length; i++) {
                y[i] = y2[i];
            }
            source2.change.emit();
        """)
    #        source2.trigger('change');



    slider2 = Slider(title="Wavelength", value=wavelist[int(len(wavelist)/2)], 
                     start=wavelist[0], end=wavelist[-1], step=sepawave,callback=callback2)


    #slider2.js_on_change('value',callback2)



    hover2 = HoverTool(
            tooltips=[
    #            ("index", "$index"),
                ("(x,y)", "($x{1}, $y)"),
            ]
        )
    c2.add_tools(hover2)
    #c2.add_tools(tools.ResizeTool())

    layout = column(slider2,c2)
    show(layout)




#Copy to other one since MEF are weird in Gemini

#sourcenameorg = 'tObj7'
#databasegeneral = 'database/apcx25sexm'
#filename = 'database/ap'+sourcenameorg
#fitsfilename = sourcenameorg+'.fits'


#sourcename = 'copy'+sourcenameorg
def getspectra(sourcenameorg,center,low,high,b_lowr=[90,70],b_upr = [143,152],xrange=(6560,6570),yrange=(0,500)):  
    
    
    center = center
    low = low
    high = high
    b_low = [center-b_lowr[0],center-b_lowr[1]]
    b_up = [b_upr[0]-center,b_upr[1]-center]
    
    #b_low = [center-90,center-70]
    #b_up = [143-center,152-center]

    
    #Databse dowlaoded from old agithub vimos cx25
    databasegeneral = 'database/apcx25sexm'
    #filename = 'database/ap'+sourcenameorg
    sourcename = 'copy'+sourcenameorg
    filename = 'database/ap'+sourcename
    #print(sourcename)
    if not os.path.isfile(sourcename+'.fits'):
        print('copy')
        iraf.imcopy(sourcenameorg+'[SCI]',sourcename)
    else:
        print('no')
    copyfile(databasegeneral,filename)
    #Now read and replace center and upper and lower values
    with open(filename) as f:
        for lines in f:
            if 'image' in lines:
                    imagenameor = lines
                    imagename= lines.replace(lines.split()[1],sourcename)
            if 'center' in lines:
                    numerocenteror = lines
                    numerocenter = lines.replace(lines.split()[2], str(center))
            if 'low' in lines:
                    numerolowor = lines
                    numerolow = lines.replace(lines.split()[2],str(low))
            if 'high' in lines:
                    numerohighor = lines
                    numerohigh = lines.replace(lines.split()[2], str(high))
            if 'xmin' in lines:
                bloworiginal  = lines
                blow = lines.replace(lines.split()[1],str(min(b_low)))
            if 'xmax' in lines:
                buporiginal  = lines
                bup = lines.replace(lines.split()[1],str(max(b_up)))
            if 'sample' in lines:
                sampleor = lines
                sampleb = '\t\tsample '+str(min(b_low))+':'+str(max(b_low))+','+str(min(b_up))+':'+str(max(b_up))+'\n'
                break

    with open(filename) as f:
        filedata = f.read()

    filedata = filedata.replace(imagenameor,imagename)
    filedata = filedata.replace(numerocenteror,numerocenter)
    filedata = filedata.replace(numerolowor,numerolow)
    filedata = filedata.replace(numerohighor,numerohigh)
    filedata = filedata.replace(bloworiginal,blow)
    filedata = filedata.replace(buporiginal,bup)
    filedata = filedata.replace(sampleor,sampleb)


    with open(filename,'w') as f:
        f.write(filedata)

    if os.path.exists(sourcename+'.ms.fits'):
        os.remove(sourcename+'.ms.fits')
    #Call them 
    iraf.noao.twodspec()
    iraf.noao.twodspec.apextract()
    iraf.noao.twodspec.apextract.setParam('dispaxis','1')
    #http://vivaldi.ll.iac.es/sieinvens/siepedia/pmwiki.php?n=HOWTOs.PythonianIRAF
    iraf.noao.apextract.apall.setParam('input',sourcename+'.fits')
    iraf.noao.apextract.apall.setParam('output',sourcename+'.ms.fits')
    iraf.noao.twodspec.apextract.apall.setParam('recenter','no')
    iraf.noao.twodspec.apextract.apall.setParam('resize','no')
    iraf.noao.twodspec.apextract.apall.setParam('edit','no')
    iraf.noao.twodspec.apextract.apall.setParam('trace','no')
    iraf.noao.twodspec.apextract.apall.setParam('interactive','no')
    iraf.noao.twodspec.apextract.apall.setParam('apertures','1')
    iraf.noao.twodspec.apextract.apall.setParam('find','no')
    iraf.noao.twodspec.apextract.apall.setParam('clean','yes')
    iraf.noao.twodspec.apextract.apall.setParam('background','average')
    iraf.noao.twodspec.apextract.apall.setParam('b_sample','-10,0:0,0')
    iraf.noao.apextract.apall.saveParList(filename='uparm/'+sourcename+'.par')
    iraf.noao.twodspec.apextract.apall(ParList='uparm/'+sourcename+'.par')
    #print(sourcename)
    
        #Plotting
    ##For srfm[0].header["CTYPE1"] = 'LINEAR'
    #Other way
    #srfm = fits.open(sourcename+'.ms.fits')
    #secondstar = srfm[0].data#[0][0]##[0][0] if using clean
    #secondstar = srfm[0].data[0][0]
    #xn = srfm[0].header["NAXIS1"]
    #refx = srfm[0].header["CRVAL1"]
    #step = srfm[0].header['CD1_1']
    #cr = srfm[0].header['CRPIX1']
    #
    #xlist = [ refx + step*(i - cr) for i in np.arange(1, len(secondstar)+1) ]
    #Create ColumnDataSource
    #x = np.array(xlist)
    #y = np.array(secondstar)
    
    
            
    ##added this becuase it is onedspec
    #if os.path.exists()
    os.rename(sourcename+'.ms.0001.fits',sourcename+'.ms.fits')


    #Not sure if I need to do this
    if os.path.exists(sourcename+'.dispcor.fits'):
        os.remove(sourcename+'.dispcor.fits')

    iraf.dispcor(sourcename+'.ms.fits',sourcename+'.dispcor.fits')
    iraf.wspectext(sourcename+'.dispcor.fits[*,1,1]',sourcename+'.txt',header='no')

    x=[]
    y=[]
    with open(sourcename+'.txt') as f:
        for lines in f:
            x.append(float(lines.split()[0]))
            y.append(float(lines.split()[1]))


    source = ColumnDataSource(data=dict(x=x,y=y))
    hover = HoverTool(
            tooltips=[
                #("index", "$index"),
                ("(x,y)", "($x{1.11}, $y)"),
            ]
        )

    #yr = (0,500)
    #xr = (6560,6570)
    xr= xrange
    yr = yrange
    plot = figure(x_axis_label='Angstrom', y_axis_label='Y',title="Spectra",
                  active_scroll='wheel_zoom',plot_width=900, plot_height=700,y_range=yr,x_range=xr)
    plot.add_tools(hover)
    #plot.add_tools(tools.ResizeTool())
    plot.line('x','y',source=source)
    show(plot)
    return x,y



def overplotlines(linesdic,x,y,xrange,yrange):
    """Put dic of lines, then the spectra to plot and the range"""


    source = ColumnDataSource(data=dict(x=x,y=y))
    hover = HoverTool(
            tooltips=[
                #("index", "$index"),
                ("(x,y)", "($x{1.11}, $y)"),
            ]
        )

    #yr = (0,500)
    #xr = (6560,6570)
    xr= xrange
    yr = yrange



    plot = figure(x_axis_label='Angstrom', y_axis_label='Y',title="Spectra",
                  active_scroll='wheel_zoom',plot_width=900, plot_height=700,y_range=yr,x_range=xr)
    plot.add_tools(hover)
    #plot.add_tools(tools.ResizeTool())
    plot.line('x','y',source=source)


    for name, xloc in linesdic.iteritems():
        yloc = y[np.where(abs(xloc - np.array(x)) < 2)[0][0]]
        span = Arrow(end=NormalHead(fill_color='orange', size=10),
                 x_start=xloc, y_start = yloc - .15, x_end = xloc, y_end= yloc-.03
                )
        plot.add_layout(span)
        my_label = Label(x=xloc, y=yloc-.03, text=name.name)
        plot.add_layout(my_label)

    show(plot)

#Parametest for the fit:
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]

def isDigit(x):
    try:
        float(x)
        return True
    except ValueError:
        return False
    
def gaussian(x, mu, sig,core):
    return core*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))





def plotgaussian(sourcename,x,y,xrange,yrange,linestofit,regionlines):
    """List linestofit = []"""

    s = sourcename

    x = np.array(x)


    line = linestofit
    regionfit = regionlines



    #Initialize files log and lines
    #! echo '' > fited.log
    open('fited.log', 'w').close()



    #Plot
    #sourceg = ColumnDataSource(data=dict(x=wavex,y=yg))


    source = ColumnDataSource(data=dict(x=x,y=y))
    hover = HoverTool(
            tooltips=[
                #("index", "$index"),
                ("(x,y)", "($x{1.11}, $y)"),
            ]
        )

    #yr = (0,500)
    #xr = (6560,6570)
    xr= xrange
    yr = yrange
    plot = figure(x_axis_label='Angstrom', y_axis_label='Y',title="Spectra",
                  active_scroll='wheel_zoom',plot_width=900, plot_height=700,y_range=yr,x_range=xr)
    plot.add_tools(hover)
    #plot.add_tools(tools.ResizeTool())
    plot.line('x','y',source=source)



    for lines in zip(line,regionfit):
        regionf = "{} {}".format(lines[1][0],lines[1][1])
        #wavelenght
        xlimns = [find_nearest(x,i)[0] for i in lines[1]    ]
        wavex = x[xlimns[0]:xlimns[1]]
        lineszero = lines[0]
        print(lineszero)
        #! echo '$lineszero' > lines.lines
        with open('lines.lines','w') as file:
            #print(lineszero)
            file.write(lineszero)



        iraf.fitprofs('copy'+s+'.ms.fits',pos='lines.lines', reg=regionf ,
                      fitbackground= 'yes', 
                      logfile='fited.log')
                      #,nerrsample='100',sigma0='4',invgain='4')


        #Plotting the gaussian
        #Find in log file  
        npattern = re.compile('[-\d.]+')
        gparameters=[]
        with open('fited.log','r') as file:
            for lines in file:
                if '(' not in lines:
                    temp = npattern.findall(lines)
                if len(temp) == 7 and all(isDigit(i) for i in temp):
                    gparameters.append(temp)

        #gaussian
        gparamfinal = [ float(i) for i in gparameters[-1] ]
        centerg, contg, fluxg, eqwg, coreg, fwhmg, fwhml = gparamfinal
        yg = gaussian(wavex,centerg,fwhmg/2.3538,coreg) + contg



        plot.line(wavex,yg,color='red')

    show(plot)

def plotspectra(fitsfilename,xrange,yrange):
        #Not sure if I need to do this
        
    filetocreate = 'e'+fitsfilename+'.fits'
    if os.path.exists(filetocreate):
        os.remove(filetocreate)

    iraf.dispcor(fitsfilename+'.ms.fits',filetocreate)
    iraf.wspectext(filetocreate+'[*,1,1]',fitsfilename+'.txt',header='no')

    x=[]
    y=[]
    with open(fitsfilename+'.txt') as f:
        for lines in f:
            x.append(float(lines.split()[0]))
            y.append(float(lines.split()[1]))


    source = ColumnDataSource(data=dict(x=x,y=y))
    hover = HoverTool(
            tooltips=[
                #("index", "$index"),
                ("(x,y)", "($x{1.11}, $y)"),
            ]
        )

    #yr = (0,500)
    #xr = (6560,6570)
    xr= xrange
    yr = yrange
    plot = figure(x_axis_label='Angstrom', y_axis_label='Y',title="Spectra",
                  active_scroll='wheel_zoom',plot_width=900, plot_height=700,y_range=yr,x_range=xr)
    plot.add_tools(hover)
    #plot.add_tools(tools.ResizeTool())
    plot.line('x','y',source=source)
    show(plot)
    #return x,y