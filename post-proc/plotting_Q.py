#!/env/python
#
#
# script plotting parameter uncertainty bands in discharge 
#
# Matthias Zink, Nov 2012
#
import numpy       as np    # module for number crunching
import scipy.stats.mstats as scim
import ufz
import matplotlib as mpl
from matplotlib.patches import Rectangle
#
######################
# input dirs and files
######################
#
QFile     = '../old_code/sub_00020/output_b1/daily_discharge.out' # 'daily_Q.out'

# -------------------------------------------------------------------------
# Command line arguments
#
pdffile   = ''
import optparse
parser = optparse.OptionParser(usage='%prog [options]',
                               description="Plotting file following template of MC.")
parser.add_option('-p', '--pdffile', action='store', dest='pdffile', type='string',
                  default=pdffile, metavar='File',
                  help='Name of pdf output file (default: open X-window).')
(opts, args) = parser.parse_args()

pdffile  = opts.pdffile
del parser, opts, args
#if os.path.exists(pdffile):   # delete existing outputfile before start new
#        sys.exit('ERROR: Output file does already exists!')
# -------------------------------------------------------------------------
# Customize plots
#
if (pdffile == ''):
    outtype = 'x'
else:
    outtype = 'pdf'

# Plot - paper_plots, but also all if not otherwise defined
nrow       = 5           # # of rows per figure
ncol       = 2           # # of columns per figure
hspace     = 0.12        # x-space between plots
wspace     = 0.02        # y-space between plots
textsize   = 18          # Standard text size
dt         = 4           # # of hours between tick marks on plots
dxabc      = 0.90        # % shift from left y-axis of a,b,c,... labels
dyabc      = 0.90        # % shift from lower x-axis of a,b,c,... labels
dyabcdown  = 0.05        # y-shift if abc in lower right corner
lwidth     = 1.5         # linewidth
elwidth    = 1.0         # errorbar line width
alwidth    = 1.0         # axis line width
msize      = 5.0         # marker size
mwidth     = 1.5         # marker edge width
# color: 'b'|'g'|'r'|'c'|'m'|'y'|'k'|'w'
#        'blue'|'green'|'red'|'cyan'|'magenta'|'yellow'|'black'|'white'
#        hex string '#eeefff' | RGB tuple (1,0.5,1) | html names 'burlywod', 'chartreuse', ...
#        grayscale intensity, e.g. '0.7', 'k'='0.0'
mcol1      = '#66c2a5'   # primary marker colour
mcol2      = '#fc8d62'       # color of second markers
mcol3      = '#7580b3'       # color of third markers
lcol1      = '0.0'       # primary line colour
lcol2      = '0.0'       # color of second lines
lcol3      = '0.0'       # color of third lines

llxbbox    = 0       # x-anchor legend bounding box
llybbox    = 0.2        # y-anchor legend bounding box
llrspace   = 0.02        # spacing between rows in legend
llcspace   = 0.2         # spacing between columns in legend
llhtextpad = 0.4         # the pad between the legend handle and text
llhlength  = 1.5         # the length of the legend handles
frameon    = True        # if True, draw a frame around the legend. If None, use rc
llxbbox2   = 0.60        # Tight bounding of symbol and text (w/o lines)
llhtextpad2= 0.          #                   "
llhlength2 = 1.0         #                   "
#
if (outtype == 'pdf'):
    mpl.use('PDF') # set directly after import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    # Customize: http://matplotlib.sourceforge.net/users/customizing.html
    mpl.rc('ps', papersize='a4', usedistiller='xpdf') # ps2pdf
else:
    import matplotlib.pyplot as plt
#
#mpl.rc('figure',     figsize=( 8.27, 11.69)) # a4 portrait
mpl.rc('figure',     figsize=(11.69,  8.27)) # a4 landscape
#mpl.rc('figure',     figsize=( 8, 5)) # for AGU 2012 Poster
mpl.rc('font',       **{'family':'sans-serif','sans-serif':['Helvetica']})
mpl.rc('font',       size=textsize)
mpl.rc('lines',      linewidth=lwidth, color='black')
mpl.rc('axes',       linewidth=alwidth, labelcolor='black')
mpl.rc('legend',     fontsize=textsize)
mpl.rc('path',       simplify=False) # do not remove
mpl.rc('text',       usetex=True)
mpl.rc('text.latex', unicode=True)

##############################################################################################
if (outtype == 'pdf'):
    print 'Plot PDF ', pdffile
    pdf_pages = PdfPages(pdffile)
else:
    print 'Plot X'
figsize = mpl.rcParams['figure.figsize']
ifig = 0
#
######################
# read discharge file
######################
#
QHead       = np.array(ufz.fread(QFile, header=True, skip=1, squeeze=False))
Qdata       = np.array(ufz.fread(QFile, skip=1, squeeze=False))
#
######################
# prepare data for plot
######################

# time  = Qdata[:,np.where(QHead == 't')[0]]
day    = Qdata[:,np.where(QHead == 'Day')[0]]
month = Qdata[:,np.where(QHead == 'Mon')[0]]
year  = Qdata[:,np.where(QHead == 'Year')[0]]
time  = ufz.date2dec(dy=day,mo=month,yr=year) - ufz.date2dec(dy=1,mo=1,yr=year[0])
Qobs  = Qdata[:,np.where(QHead == 'Qobs_0000411')[0]]
Qcal  = Qdata[:,np.where(QHead == 'Qsim_0000411')[0]]
#
###################################################################
# plot
###################################################################
ifig += 1
print 'Plot - Fig ', ifig
fig        = plt.figure(ifig)
nPlot = 1

ax        = fig.add_subplot(111)
#
ax.set_xlim(1, Qcal.shape[0])
ax.set_xlabel('Time')
ax.set_ylabel('Discharge Q [m$^3~s^{-1}$]')
# simulation
l1 = ax.plot(time, Qcal, 'k-', label=ur'Q$_{sim}$')
# observed series
l2 = ax.plot(time, Qobs, marker='o', linestyle='', markerfacecolor='None', markeredgecolor='r' , markeredgewidth=0.8, markersize=5, label=ur'Q$_{obs}$')
#
plt.title('First run mHM v5.0: 20.02.2013 19:12:50')
NSE = 1. - np.sum((Qcal-Qobs)**2) / np.sum((Qobs-np.mean(Qobs))**2)
plt.figtext(0.15, 0.85, r'NSE = '+ str(np.round(NSE,4)) )
# legend
handle, label = ax.get_legend_handles_labels()
ax.legend(handle, label)

#
if (outtype == 'pdf'):
    pdf_pages.savefig(fig)
    plt.close()
else:
    plt.show()

if (outtype == 'pdf'):
    pdf_pages.close()

