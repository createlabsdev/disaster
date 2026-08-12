#!/usr/bin/python3

#############################################################################
#
# MODULE:       r.scarp
#
# AUTHOR:       Martin Mergili
#
# PURPOSE:      Script for deriving landslide scarps
#               from multiple sliding planes
#
# VERSION:      20240614 (14 June 2024)
#
# COPYRIGHT:    (c) 2024 by the author
#               (c) 2024 by the University of Graz
#               (c) 2024 by the Austrian Torrent and Avalanche Control
#
#               This program is free software under the GNU General Public
#               License (>=v2). Read the file COPYING that comes with GRASS
#               for details.
#
#############################################################################

#%module
#% description: Script for deriving landslide scarps from multiple sliding planes
#% keywords: Raster
#% keywords: Topographic analysis
#% keywords: Landslide
#%end

#%option
#% key: cellsize
#% type: string
#% description: Cell size (input, m)
#% required: no
#% multiple: no
#%end

#%option
#% key: elevation
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of elevation raster map (input, m)
#% required: yes
#% multiple: no
#%end

#%option
#% key: perimeter
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of maximum scarp extent raster map (input, binary)
#% required: yes
#% multiple: no
#%end

#%option
#% key: elevscarp
#% type: string
#% gisprompt: new,raster,dcell
#% description: Name of scarp elevation raster map (output)
#% required: yes
#% multiple: no
#%end

#%option
#% key: hrelease
#% type: string
#% gisprompt: new,raster,dcell
#% description: Name of release height raster map (output)
#% required: yes
#% multiple: no
#%end

#%option
#% key: planes
#% type: string
#% description: Comma-separated data for each sliding plane (xseed, yseed, dip, slope)
#% required: no
#% multiple: yes
#%end

#Importing libraries

import grass.script as grass
from grass.script import core as grasscore
from grass.pygrass.raster import RasterSegment
import math
import sys

def main(): #starting main function

    #Reading and preparing parameters

    cellsize=options['cellsize']
    elevation=options['elevation']
    perimeter=options['perimeter']
    elevation=options['elevation']
    elevscarp=options['elevscarp']
    hrelease=options['hrelease']
    planes=options['planes']

    #Setting and reading region

    grass.run_command('g.region', flags='d')
    if cellsize: grass.run_command('g.region', flags='a', res=cellsize)

    c = grass.region()
    m = int(c['rows'])
    n = int(c['cols'])
    res = float(c['nsres'])
    north = float(c['n'])
    west = float(c['w'])

    #Preparing data

    planes=list(map(str, planes.split(',')))
    nplanes=int(len(planes)/5)

    xseed=[float(planes[0])]
    yseed=[float(planes[1])]
    zseed=[float(planes[2])]
    dip=[float(planes[3])]
    slope=[float(planes[4])]

    betax = [math.tan( slope[0] * math.pi / 180 ) * math.cos( dip[0] * math.pi / 180 )]
    betay = [math.tan( slope[0] * math.pi / 180 ) * math.sin( dip[0] * math.pi / 180 )]
   
    for i in range(1, nplanes):
    
        xseed.append(float(planes[5*i]))
        yseed.append(float(planes[5*i+1]))
        zseed.append(float(planes[5*i+2]))
        dip.append(float(planes[5*i+3]))
        slope.append(float(planes[5*i+4]))
        
        betax.append(math.tan( slope[i] * math.pi / 180 ) * math.cos( dip[i] * math.pi / 180 ))
        betay.append(math.tan( slope[i] * math.pi / 180 ) * math.sin( dip[i] * math.pi / 180 ))

    #Opening raster maps

    elev = RasterSegment(elevation)
    elev.open('r')

    peri = RasterSegment(perimeter)
    peri.open('r')

    scarp = RasterSegment(elevscarp)
    scarp.open('w', 'DCELL', overwrite=True)

    hrel = RasterSegment(hrelease)
    hrel.open('w', 'DCELL', overwrite=True)

    #Computing scarp and release height

    for i in range(0, m):
        for j in range(0, n):
        
            scarp[i, j] = zseed[0] + ( north - float(i) * res - yseed[0] ) * betay[0] - ( west + float(j) * res - xseed[0] ) * betax[0]

            for k in range(1, nplanes):
 
                 scarp[i, j] = max( scarp[i][j], zseed[k] + ( north - float(i) * res - yseed[k] ) * betay[k] - ( west + float(j) * res - xseed[k] ) * betax[k] )

            if peri[i, j] > 0: 
            
                scarp[i, j] = min( elev[i][j], scarp[i][j] )
                
            else:
                
                scarp[i, j] = elev[i][j]

    #Writing output and closing raster maps

    for i in range(0, m):
        for j in range(0, n):

            hrel[i, j] = elev[i][j] - scarp[i][j]

    elev.close()
    peri.close()
    scarp.close()
    hrel.close()

    grass.run_command('g.region', flags='d')

    sys.exit()

if __name__=='__main__':
    options, flags=grass.parser()
    main()
