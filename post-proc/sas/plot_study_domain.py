#!/usr/bin/env python
#
# purpose: creating netdcf for small subbasins in catch having the same
# georeferencing (lower left corner)
#
# input: catch mask, lut subbasin masks, subbasin mask (ASCII files)
#
# restrictions: only integer values in all input data excpected/accepted
#
# created by Matthias Zink Feb. 2014
#
import numpy as np                       # array manipulation
import ufz
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from cartopy.feature import ShapelyFeature
from cartopy.io import srtm
#
# -------------------------------------------------------------------------
# Command line arguments
#
pdffile  = ''
pngfile  = ''
import optparse
parser  = optparse.OptionParser(usage='%prog [options]',
                               description="Plotting file following template of MC.")
parser.add_option('-p', '--pdffile', action='store', dest='pdffile', type='string',
                  default=pdffile, metavar='File',
                  help='Name of pdf output file (default: open X-window).')
parser.add_option('-g', '--pngfile', action='store',
                    default=pngfile, dest='pngfile', metavar='pngfile',
                    help='Name basis for png output files (default: open screen window).')
(opts, args) = parser.parse_args()

pdffile  = opts.pdffile
pngfile  = opts.pngfile
del parser, opts, args

if (pdffile != '') & (pngfile != ''):
    raise ValueError('PDF and PNG are mutually exclusive. Only either -p or -g possible.')

# -------------------------------------------------------------------------
# Customize plots
#
if (pdffile == ''):
    if (pngfile == ''):
        outtype = 'x'
    else:
        outtype = 'png'
else:
    outtype = 'pdf'


# Plot - paper_plots, but also all if not otherwise defined
nrow       = 1           # # of rows per figure
ncol       = 1           # # of columns per figure
hspace     = 0.01        # x-space between plots
wspace     = 0.01        # y-space between plots
figleft    = 0.14        # left border ogf the figure
figright   = 0.97        # right border ogf the figure
figbottom  = 0.05        # lower border ogf the figure
figtop     = 0.95        # upper border ogf the figure
textsize   = 18          # Standard text size
dt         = 4           # # of hours between tick marks on plots
dxabc      = 0.90        # % shift from left y-axis of a,b,c,... labels
dyabc      = 0.90        # % shift from lower x-axis of a,b,c,... labels
dyabcdown  = 0.05        # y-shift if abc in lower right corner
lwidth     = 0.5         # linewidth
elwidth    = 1.0         # errorbar line width
alwidth    = 1.0         # axis line width
msize      = 10.         # marker size
mwidth     = 4           # marker edge width
bxwidth    = 0.85        # boxlplot width
# color: 'b'|'g'|'r'|'c'|'m'|'y'|'k'|'w'
#        'blue'|'green'|'red'|'cyan'|'magenta'|'yellow'|'black'|'white'
#        hex string '#eeefff' | RGB tuple (1,0.5,1) | html names 'burlywod', 'chartreuse', ...
#        grayscale intensity, e.g. '0.7', 'k'='0.0'
mcol1      = '#67A9CF'   # color of second markers
mcol2      = '#A1D99B'   # color of third markers
mcol3      = '#EF8A62'         # primary line colour
mcol4      = 'r'
lcol2      = '0.5'       # color of second lines
lcol3      = '0.0'       # color of third lines

llxbbox    = 0.5        # x-anchor legend bounding box
llybbox    = 0.87        # y-anchor legend bounding box
llrspace   = 0.02        # spacing between rows in legend
llcspace   = 1.0         # spacing between columns in legend
llhtextpad = 0.2         # the pad between the legend handle and text
llhlength  = 0.9         # the length of the legend handles
frameon    = True        # if True, draw a frame around the legend. If None, use rc
llxbbox2   = 0.60        # Tight bounding of symbol and text (w/o lines)
llhtextpad2= 0.          #                   "
llhlength2 = 1.0         #                   "

# PNG
dpi         = 300
transparent = False
bbox_inches = 'tight'
pad_inches  = 0.02
#
if (outtype == 'pdf'):
    mpl.use('PDF') # set directly after import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    # Customize: http://matplotlib.sourceforge.net/users/customizing.html
    mpl.rc('ps', papersize='a4', usedistiller='xpdf') # ps2pdf
elif (outtype == 'png'):
    mpl.use('Agg') # set directly after import matplotlib
    import matplotlib.pyplot as plt
    mpl.rc('text.latex', unicode=True)
    mpl.rc('savefig', dpi=dpi, format='png')
else:
    import matplotlib.pyplot as plt
#
# scale figsize for pdf
latexwidth    = 470                   # invoke latex document textwidth with \the\textwidth
factor        = 1.00                  # factor used in latex document like width=factor\textwidth

fig_width_in  = latexwidth * factor * 0.01388888  # figure width in inches * inches per point
fig_height_in = fig_width_in #* 1.318   # figure height in inches * golden ratio

mpl.rc('figure', figsize=( fig_width_in, fig_height_in)) # half side A4potrait, golden ratio

#mpl.rc('figure',     figsize=( 8, 5)) # for AGU 2012 Poster
#for paper half of the size because of 2 columns
# one column of 2 column paper
#mpl.rc('figure',     figsize=( 8.27/2., 8.27/2./1.618)) # half side A4potrait, golden ratio
#mpl.rc('figure',     figsize=( 8.27, 11.69)) # a4 portrait
#mpl.rc('figure',     figsize=(11.69,  8.27)) # a4 landscape
#mpl.rc('figure',     figsize=( 8.27, 11.69/2)) # a4 portrait
mpl.rc('font',       **{'family':'serif','serif':['times']})
mpl.rc('font',       size=textsize)
mpl.rc('legend',     fontsize=textsize)
mpl.rc('lines',      linewidth=lwidth, color='black')
mpl.rc('axes',       linewidth=alwidth, labelcolor='black')
mpl.rc('path',       simplify=False) # do not remove
mpl.rc('text',       usetex=True)
mpl.rc('text.latex', unicode=True)

##############################################################################################
if (outtype == 'pdf'):
    print 'Plot PDF ', pdffile
    pdf_pages = PdfPages(pdffile)
elif (outtype == 'png'):
    print('Plot PNG ', pngfile)
else:
    print 'Plot X'

figsize = mpl.rcParams['figure.figsize']
ifig = 0

shpf_eu_adm        = '../wb_paper_de/study_domain/data/eu_admin_bound.shp'
shpf_unstr         = '../wb_paper_de/study_domain/data/BuLaender.shp'
shpf_river         = '../wb_paper_de/study_domain/data/neckar_mulde_saale.shp'
shpf_river_unstrut = './shape/unstrut_river.shp'
shpf_basin_unstrut = './shape/unstrut.shp'
shpf_thubasin      = './shape/watershed_WGS84.shp'
shpf_thubasin      = './shape/Grenzen_WGS84'
txt_ecstat         = './EC_wb_paper.txt'

demfile            = '../wb_paper_de/study_domain/data/dem_de/dem_1km.nc'

inpath             = '../wb_paper_de/study_domain/data/watersheds/'

catchlist   = ['saale','weser']
label_cat   = ['Saale','Weser']
lon_cat     = [  10.90,   9.20]
lat_cat     = [  51.85,  52.40]
# MAPPLOT #########################################
ifig          += 1
print 'Plot - Fig ', ifig
fig            = plt.figure(ifig)
iplot          = 1

# data PROCESSING #############################################################
data          = ufz.readnetcdf(demfile, var='mask')[0,:,:]
lons          = ufz.readnetcdf(demfile, var='lon')
lats          = ufz.readnetcdf(demfile, var='lat')
dem           = np.ma.array(data, mask=~(data>-9999.))
shaded        = np.ma.array(srtm.add_shading(data*10,315.,45.), mask=~(data>-9999.))
test          = dem / 5. + shaded
smooth        = np.ma.array(ufz.savitzky_golay2d(shaded,7,2), mask=~(data>-9999.))

# ec stations
data        = ufz.sread(txt_ecstat, skip=1, strarr=True)
data_header = ufz.sread(txt_ecstat, skip=1, header=True, strarr=True)

lat_ec      = data[:,np.where(data_header=='Latitude')[0][0]].astype('float')
lon_ec      = data[:,np.where(data_header=='Longitude')[0][0]].astype('float')
label_ec    = data[:,np.where(data_header=='Stat_ID')[0][0]]
id_ec       = data[:,np.where(data_header=='Id')[0][0]]
# plotting    #############################################################

mp      = fig.add_axes(ufz.position(nrow,ncol,1, wspace=0.01, hspace=0.01, 
                                    left=figleft, right=figright, top=figtop, bottom=figbottom), 
                       projection=ccrs.Mercator())

# coordinate transformation
transform = ccrs.PlateCarree()._as_mpl_transform(mp)
# coloring
#cmap          = 'binary'
#cmap2         = 'gist_gray'
colors        = ufz.get_brewer('gsdtol',rgb=True)[1:]
cmap          = mpl.colors.ListedColormap(colors[:4:-1])
#cmap          = mpl.colors.ListedColormap(ufz.get_brewer('OceanLakeLandSnow',rgb=True)[1:])
cmap2         = mpl.colors.ListedColormap(colors)

breaks        = list(np.linspace(-50,1000,22)) + [1200,1400,1600,1800,2000,2200]
norm          = mpl.colors.BoundaryNorm(breaks, cmap.N)
norm2         = mpl.colors.BoundaryNorm(np.linspace(-20,250,100),cmap.N)

# transparency
cmap2._init()
cmap2._lut[:-3:,-1] = 0.99
cmap._init()
cmap._lut[:-3:,-1]  = 0.015


#pcm2    = mp.pcolormesh(lons, lats, test,   transform=ccrs.PlateCarree(), cmap=cmap, norm=norm2, 
#                        edgecolors='face')
# pcm1    = mp.pcolormesh(lons, lats, smooth, transform=ccrs.PlateCarree(), cmap=cmap2,
#                         edgecolors='face')#, antialiased=True)
# pcm3    = mp.pcolormesh(lons, lats, dem,    transform=ccrs.PlateCarree(), cmap=cmap, norm=norm, 
#                         edgecolors='face') # alpha=0.08, 
#pcm3    = mp.contourf(lons, lats, dem,   transform=ccrs.PlateCarree(), 
#                      levels=breaks, cmap=cmap, norm=norm, alpha=0.5)#, antialiased=True)

#plt.show()
#stop
# catchments
for icat in range(len(catchlist)): 
    # federal countries, bundeslaender
    if   (icat==5 or icat==2):
        color='r'
    elif (icat==0 or icat==4):
        color='g'
    else:
        color='b'
    reader    = shpreader.Reader(inpath + 'basin_' + catchlist[icat] + '.shp')
    fed_count = ShapelyFeature(reader.geometries(), ccrs.PlateCarree())
    mp.add_feature(fed_count, facecolor=color, alpha=0.10, lw=0.5, edgecolor='0') #, alpha=0.75)
    # catchment labels
    mp.annotate(label_cat[icat],(lon_cat[icat], lat_cat[icat]), xycoords=transform,
                va='center', ha='left', color=color)


# eddy stations
mp.plot(lon_ec, lat_ec, ls='None',marker='s',mfc='None', mec='k',mew=1.2, ms=3,transform=ccrs.PlateCarree())
for i in range(len(lon_ec)):
    if ( label_ec[i] == 'E1'):
        xp=lon_ec[i]+0.03; yp=lat_ec[i]#-0.15
        va='top'; ha='left'
    elif ( label_ec[i] == 'E2'):
        xp=lon_ec[i]; yp=lat_ec[i]#-0.15
        va='top'; ha='center'
    elif ( label_ec[i] == 'E3'):
        xp=lon_ec[i]-0.1; yp=lat_ec[i]
        va='bottom'; ha='right'
    elif ( label_ec[i] == 'E4'):
        xp=lon_ec[i]; yp=lat_ec[i]#-0.15
        va='top'; ha='right'
    mp.annotate(id_ec[i].split('-')[1],(xp, yp), xycoords=transform, va=va, ha=ha, fontsize=0.95*textsize)
    #plt.text(lon_ec[i], lat_ec[i], label_ec[i], transform=ccrs.PlateCarree())

# mp.set_global()
ixticks  = np.arange(  2,17.5, 1)
iyticks  = np.arange( 46,  60, 1)

mp.set_yticks(iyticks,  crs=ccrs.PlateCarree()) 
mp.set_xticks(ixticks,  crs=ccrs.PlateCarree())

xnames = []
for i in list(ixticks):
  if (i < 0):
      xnames.append(r'$\mathrm{'+str(int(i))+'\,^{\circ} W}$')
  elif (i == 0):
      xnames.append(r'$\mathrm{'+str(int(i))+'\,^{\circ}}$')
  elif (i > 0):
      xnames.append(r'$\mathrm{'+str(int(i))+'\,^{\circ} E}$')
ynames     = [ r'$\mathrm{'+str(i)+'\,^{\circ} N}$' for i in iyticks]

yticknames = plt.setp(mp, yticklabels=ynames)
plt.setp(yticknames, fontsize=.9*textsize) # fontsize='medium')
xticknames = plt.setp(mp, xticklabels=xnames)
plt.setp(xticknames, fontsize=.9*textsize)

# properties and features
mp.set_extent([7.95, 12.95, 50.1, 53.95]) #, crs=ccrs.Mercator())

# geographical features
ocean         = cfeature.NaturalEarthFeature(category='physical', name='ocean',
                                             scale='50m', facecolor=cfeature.COLORS['water'])
mp.add_feature(ocean, edgecolor='none')
# state_borders = cfeature.NaturalEarthFeature(category='cultural', name='admin_0_countries',
#                                              scale='10m', facecolor='none')
# mp.add_feature(state_borders, edgecolor='k', lw=0.5)
reader    = shpreader.Reader(shpf_eu_adm)
state_borders = ShapelyFeature(reader.geometries(), ccrs.PlateCarree())
mp.add_feature(state_borders, facecolor='none', edgecolor='k', lw=0.5)
# rivers
rivers = cfeature.NaturalEarthFeature(category='physical', name='rivers_lake_centerlines',
                                      scale='10m', facecolor='none')
mp.add_feature(rivers, edgecolor='#045a8d', lw=0.5) 

reader    = shpreader.Reader(shpf_river)
add_rivers = ShapelyFeature(reader.geometries(), ccrs.PlateCarree())
mp.add_feature(add_rivers, facecolor='none', edgecolor='#045a8d', lw=0.5)

reader        = shpreader.Reader(shpf_basin_unstrut)
unstrut_basin = ShapelyFeature(reader.geometries(), ccrs.PlateCarree())
mp.add_feature(unstrut_basin, facecolor='0.5', edgecolor='0.5', lw=1.0)

reader        = shpreader.Reader(shpf_thubasin)
thube_basin  = ShapelyFeature(reader.geometries(), ccrs.PlateCarree())
mp.add_feature(thube_basin, facecolor='none', edgecolor='k', lw=1.0)

reader        = shpreader.Reader(shpf_river_unstrut)
river_unstrut = ShapelyFeature(reader.geometries(), ccrs.PlateCarree())
mp.add_feature(river_unstrut, facecolor='none', edgecolor='#045a8d', lw=0.5)


if (outtype == 'pdf'):
  pdf_pages.savefig(fig)
  plt.close()
elif (outtype == 'png'):
  fig.savefig(pngfile+str(ifig).zfill(3)+'.png', transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
  plt.close(fig)
else:
  plt.show()

if (outtype == 'pdf'):
    pdf_pages.close()
elif (outtype == 'png'):
    pass


