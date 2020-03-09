
import os
import csv
import glob
import matplotlib.pyplot as plt
from operator import itemgetter
import MoNeT_MGDrivE.plots as plots
import MoNeT_MGDrivE.spatialVideosLegacy as video
from mpl_toolkits.basemap import Basemap

PAD = .1
COLORS = [
        plots.rescaleRGBA((47, 28, 191, 255/2.5)),    # 0: Faded navy blue
        plots.rescaleRGBA((0, 169, 255, 255/7.5)),    # 1: Cyan
        plots.rescaleRGBA((255, 0, 152, 255/1)),      # 2: Magenta
        plots.rescaleRGBA((37, 216, 17, 255/6)),      # 3: Bright green
        plots.rescaleRGBA((255, 255, 255, 255/1)),    # 4: White
        plots.rescaleRGBA((0, 0, 0, 255/5))           # 5: Black
    ]


def get_corners(fileName):
    lats = []
    longs = []
    clusterData = open(fileName, 'r')
    next(clusterData)
    for line in clusterData:
        tokens = line.split(',')
        (lat, long) = (float(tokens[1]), float(tokens[2]))
        lats.append(lat)
        longs.append(long)

    minLat = min(lats)
    minLong = min(longs)
    maxLat = max(lats)
    maxLong = max(longs)
    return [[minLong, maxLong], [minLat, maxLat]]


def createBasemapInstance(minLat, maxLat, minLon, maxLon, pad=1.5):
    base = Basemap(
            projection='merc',
            lat_0=(maxLat - minLat)/2, lon_0=(maxLon - minLon)/2,
            resolution='h', area_thresh=0.1,
            llcrnrlon=minLon - pad, llcrnrlat=minLat - pad,
            urcrnrlon=maxLon + pad, urcrnrlat=maxLat + pad,
            epsg=4269
        )
    return base


def populateClustersFromList(
            cList, pFileLocation, pFilePattern={}
        ):
    # Create the empty list to contain the clusters
    (clusterNum, clusters) = (len(set(cList)), [])
    for i in range(clusterNum):
        clusters.append({'male': [], 'female': []})
    # Male files clustering
    if 'male' in pFilePattern:
        patchFileList = sorted(glob.glob(pFileLocation+pFilePattern['male']))
    else:
        patchFileList = sorted(glob.glob(pFileLocation+'/M_*'))
    for index, patchFileN in enumerate(patchFileList):
        clusters[cList[index]]['male'].append(patchFileN)
    # Female files clustering
    if 'female' in pFilePattern:
        patchFileList = sorted(glob.glob(pFileLocation+pFilePattern['female']))
    else:
        patchFileList = sorted(glob.glob(pFileLocation+'/F_*'))
    for index, patchFileN in enumerate(patchFileList):
        clusters[cList[index]]['female'].append(patchFileN)
    # Return
    return clusters


def draw_dots(m, alphas, colorList, long=0, lat=0, size=60):
    # start = 0.0
    for idx, value in enumerate(alphas):
        m.scatter(
                [long], [lat], latlon=True, marker=(6, 0),
                s=max(6, 0.11 * size), facecolor=colorList[idx],
                alpha=value, linewidths=.25, edgecolors='White'
            )


def generateClusterGraphs(
            clstFile,
            aggList, coordinates, destination, colorList, original_corners,
            padding, dpi, countries=False, skip=False, refPopSize=1,
            verbose=True, background=False, timeLocation=(.5, .5),
            colors=COLORS
        ):
    time = len(aggList[0])
    timeMax = list(range(time))
    for tick in timeMax:
        imgFileName = destination+'/c_'+str(tick).zfill(6)+".png"
        if skip and os.path.isfile(imgFileName):
            continue

        for idx, cData in enumerate(aggList):
            if idx == 0:
                (fig, ax, m) = createMap(clstFile, COLORS, pad=.025)
            pops = []
            try:
                pops = cData[tick]
                alphas, size = video.getAlphas(pops)
                if alphas:
                    draw_dots(
                            m, alphas, colorList,
                            coordinates[1][idx], coordinates[0][idx],
                            size/refPopSize
                        )
                else:
                    continue
            except Exception as e:
                return e
        else:
            ax.axis('off')
            plt.text(
                    timeLocation[0], timeLocation[1], str(tick+1).zfill(4),
                    ha='left', va='top',
                    transform=fig.transFigure
                )
            fig.savefig(imgFileName,
                        dpi=dpi, orientation='portrait', papertype=None,
                        transparent=False, format="png",
                        bbox_inches='tight', pad_inches=0.05, frameon=None)
            plt.close(fig)
            plt.close('all')
            if original_corners:
                fig, ax, m = video.createFig(
                        original_corners, padding, countries
                    )
            else:
                fig, ax, m = video.createFig(coordinates, padding, countries)
        if verbose:
            print(
                    '* Exporting frame ({}/{})'.format(
                            str(tick+1).zfill(5), str(time).zfill(5)
                    ), end='\r'
                )
    return


def createMap(clusterFile, COLORS, pad=.025):
    (minLat, maxLat, minLong, maxLong) = (0, 0, 0, 0)
    (lats, longs, clusters) = ([], [], [])
    clusterData = open(clusterFile, 'r')
    next(clusterData)
    for line in clusterData:
        tokens = line.split(',')
        (lat, long, cluster) = (
                float(tokens[1]), float(tokens[2]), int(tokens[3])
            )
        lats.append(lat)
        longs.append(long)
        clusters.append(cluster)

    (minLat, maxLat, minLong, maxLong) = (
            min(lats), max(lats), min(longs), max(longs)
        )
    (minCluster, maxCluster) = (min(clusters), max(clusters))

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111, label="1")
    m = Basemap(
            projection='merc',
            llcrnrlat=minLat-pad, urcrnrlat=maxLat+pad,
            llcrnrlon=minLong-pad, urcrnrlon=maxLong+pad,
            lat_ts=20, resolution='i', ax=ax
        )
    m.drawcoastlines(color=COLORS[1], linewidth=5, zorder=-1)
    m.drawcoastlines(color=COLORS[0], linewidth=2, zorder=-1)
    m.drawcoastlines(color=COLORS[1], linewidth=.5, zorder=-1)
    # m.fillcontinents(color=COLORS[3], lake_color='aqua')
    m.scatter(
            longs, lats, latlon=True, alpha=.1, marker='x', s=1,
            cmap=plt.get_cmap('winter'), c=clusters,
            vmin=minCluster, vmax=maxCluster
        )
    ax.tick_params(
            axis='both',       # changes apply to the both
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            left=False,
            right=False,
            labelbottom=False,  # labels along the bottom edge are off
            labelleft=False
        )
    ax.axis('off')
    return (fig, ax, m)


def getClustersFromAggFiles(coordinatesFileI):
    coordinates = []
    clusterFile = open(coordinatesFileI, 'r')
    for (i, line) in enumerate(clusterFile):
        if i > 0:
            tokens = line.split(',')
            coordinates.append((
                    float(tokens[4]),
                    float(tokens[5]),
                    int(tokens[3])
                ))
    coordinates = list(set(coordinates))
    sortedCoordinates = sorted(coordinates, key=itemgetter(2))
    (lats, lons) = (
            [i[0] for i in sortedCoordinates],
            [i[1] for i in sortedCoordinates]
        )
    return (lats, lons)


def readClustersIDs(filepath):
    with open(filepath) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        clustersList = []
        for (i, row) in enumerate(csv_reader):
            if i != 0:
                clustersList.append(int(row[3]))
    return clustersList
