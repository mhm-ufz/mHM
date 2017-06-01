#from pylab import *
import ufz
import matplotlib as mpl
import matplotlib.pyplot as plt
#from ufz import readnetcdf
from sas import *

from cartopy.feature import ShapelyFeature
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader

def plot_routing_network(lon, lat, routing_array):
    
    ax = axes()

#    ax = axes( projection=ccrs.Mercator())
#    ax.pcolormesh(lon, lat, routing_array[ :, :], transform=ccrs.PlateCarree(), cmap='RdYlGn', edgecolors='face')
#    ax.gridlines(draw_labels=True)
    
#    ax.arrow(0, 0, 1.5, 0.5, head_width=0.05, head_length=0.1, fc='k', ec='k')
    ax.annotate('axes center', xy=(.5, .5),  xycoords='axes fraction', 
                horizontalalignment='center', verticalalignment='center')
    
#    fig = figure(1)
#    ax = fig.add_axes(ufz.position(1,1,1, wspace=0.01, hspace=0.01, left=0.14, right=0.97, top=0.95, bottom=0.05), projection=ccrs.Mercator())
    

    show()
    
def plot_basin(lon, lat, basin_data, path = None, color_map = None):
    
    mpl.rc('text',       usetex=True)
    mpl.rc('text.latex', unicode=True)
    
#    shpf_river_unstrut = './shape/unstrut_river.shp'
#    shpf_basin_unstrut = './shape/unstrut.shp'
#    shpf_thubasin      = './shape/watershed_WGS84.shp'
#    shpf_thubasin      = './shape/Grenzen_WGS84'
    
#    if path:
#        x = readnetcdf(path, var='easting')
#        y = readnetcdf(path, var='northing')   
    
#        nagelstedt = MapPoint(x, y, lon, lat, 4409672, 5664533)
#        cammerbach = MapPoint(x, y, lon, lat, 4392737, 5668127)
#        suthbach = MapPoint(x, y, lon, lat, 4396321, 5667869)
        
#        reckenbuhl = MapPoint(x, y, lon, lat, 4387837, 5664518)
#        bechstedt = MapPoint(x, y, lon, lat, 4388942, 5664961)
#        heuberg = MapPoint(x, y, lon, lat, 4390215, 5665475)
#        fuchsloch = MapPoint(x, y, lon, lat, 4391430, 5665590)
#        ruspelsweg = MapPoint(x, y, lon, lat, 4392925, 5666035)
    
#    test = np.append(lon.data, np.zeros((1, 90)), 0) 
#    print(test.shape)

    fig = plt.figure(1)

    font = {'size' : 36}
    mpl.rc('font', **font)
#    cmap = mpl.cm.get_cmap('YlOrRd')
    mp = fig.add_axes(ufz.position(1,1,1), projection=ccrs.Mercator())
    pcm = mp.pcolormesh(lon, lat, basin_data, transform=ccrs.PlateCarree(), 
                        cmap=color_map, edgecolors='face')
    ixticks  = np.round(np.arange((lon.min()), (lon.max()), 0.1)*10)/10
    iyticks  = np.round(np.arange((lat.min()), (lat.max()), 0.1)*10)/10

    mp.set_yticks(iyticks,  crs=ccrs.PlateCarree()) 
    ynames     = [ r'$\mathrm{'+str((i))+'\,^{\circ} N}$' for i in iyticks]
    yticknames = plt.setp(mp, yticklabels=ynames)
    mp.set_xticks(ixticks,  crs=ccrs.PlateCarree())
    xnames = []
    for i in list(ixticks):
        if   (i < 0):
            xnames.append(r'$\mathrm{'+str((i))+'\,^{\circ} W}$')
        elif (i == 0):
            xnames.append(r'$\mathrm{'+str((i))+'\,^{\circ}}$')
        elif (i > 0):
            xnames.append(r'$\mathrm{'+str((i))+'\,^{\circ} E}$')
    xticknames = plt.setp(mp, xticklabels=xnames)
#    plt.axis([lon.min(), lon.max(), lat.min(), lat.max()])
#    plt.xticks(ixticks)
#    mp.tick_params(labelbottom=True, labelleft=False)  
    mp.grid(True)                  
    cb = plt.colorbar(pcm, orientation='vertical')
#    mpl.rcParams['xtick.labelsize'] = font_size 
    
#    fig = plt.figure(1)
#    ax = fig.add_axes(ufz.position(1,1,1, wspace=0.01, hspace=0.01, left=0.14, right=0.97, top=0.95, bottom=0.05), projection=ccrs.Mercator())

#    reader = shpreader.Reader(shpf_river_unstrut)
#    riv_unstrut = ShapelyFeature(reader.geometries(), ccrs.PlateCarree())
#    ax.add_feature(riv_unstrut, facecolor='none', edgecolor='#045a8d', lw=2.5)

#    reader = shpreader.Reader(shpf_basin_unstrut)
#    unstrut_basin = ShapelyFeature(reader.geometries(), ccrs.PlateCarree())
#    ax.add_feature(unstrut_basin, facecolor='none', edgecolor='0.5', lw=1.0)
    
#    reader = shpreader.Reader(shpf_thubasin)
#    thube_basin = ShapelyFeature(reader.geometries(), ccrs.PlateCarree())
#    ax.add_feature(thube_basin, facecolor='none', edgecolor='k', lw=1.0)
    
#    if path:
##        ax.scatter(nagelstedt.lon, nagelstedt.lat, transform=ccrs.PlateCarree())
#        ax.scatter(reckenbuhl.lon, reckenbuhl.lat, c='b', 
#                   transform=ccrs.PlateCarree())
#        ax.scatter(bechstedt.lon,  bechstedt.lat,  c='g', 
#                   transform=ccrs.PlateCarree())
#        ax.scatter(heuberg.lon,    heuberg.lat,    c='r', 
#                   transform=ccrs.PlateCarree())
#        ax.scatter(fuchsloch.lon,  fuchsloch.lat,  c='c', 
#                   transform=ccrs.PlateCarree())
#        ax.scatter(ruspelsweg.lon, ruspelsweg.lat, c='m', 
#                   transform=ccrs.PlateCarree())
    
#    ax.set_autoscaley_on(False)
#    ax.set_ylim([0,10])
#    ax.set_extent( [ 9.5, 12.5, 50.5, 51.7 ] )
    
    
#    ax.title('mean travel time', fontsize=font_size)
    
    plt.show()    